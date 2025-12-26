#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>
using namespace std;

const double pixelareacon = 0.09;
const double weightcarbon = 0.35;
const double weightsloper = 0.25;
const double weightproxim = 0.25;
const double weightedgext = 0.15;

vector<float> normalizgrid(const vector<float>& inputgridval,
                           float minboundvalu,
                           float maxboundvalu) {
    vector<float> normgridvalu(inputgridval.size());

    for(size_t gridindexval = 0; gridindexval < inputgridval.size(); ++gridindexval) {
        float currentvalue = inputgridval[gridindexval];

        if(currentvalue < minboundvalu) currentvalue = minboundvalu;
        if(currentvalue > maxboundvalu) currentvalue = maxboundvalu;

        normgridvalu[gridindexval] =
            (currentvalue - minboundvalu) / (maxboundvalu - minboundvalu);
    }
    return normgridvalu;
}

static inline void disttransone(const float* sourcevalues,
                                int valuecount01,
                                float* distancevalu,
                                int* vertexindexs,
                                float* breakpointev) {
    const float infinitevalu = numeric_limits<float>::infinity();

    int currentindex = 0;
    vertexindexs[0]  = 0;
    breakpointev[0]  = -infinitevalu;
    breakpointev[1]  =  infinitevalu;

    for(int queryindexer = 1; queryindexer < valuecount01; ++queryindexer) {
        float intersectpos;

        while (true) {
            const int vertexpicked = vertexindexs[currentindex];

            intersectpos =
                ((sourcevalues[queryindexer] +
                  queryindexer * (float)queryindexer) -
                 (sourcevalues[vertexpicked] +
                  vertexpicked * (float)vertexpicked)) /
                (2.0f * (queryindexer - vertexpicked));

            if(intersectpos > breakpointev[currentindex]) {
                break;
            }

            if(--currentindex < 0) {
                currentindex = 0;
                break;
            }
        }

        ++currentindex;
        vertexindexs[currentindex]   = queryindexer;
        breakpointev[currentindex]   = intersectpos;
        breakpointev[currentindex+1] = infinitevalu;
    }

    currentindex = 0;

    for(int queryindexer = 0; queryindexer < valuecount01; ++queryindexer) {
        while (breakpointev[currentindex + 1] < queryindexer) {
            ++currentindex;
        }

        float offsetvaluse =
            queryindexer - (float)vertexindexs[currentindex];

        distancevalu[queryindexer] =
            offsetvaluse * offsetvaluse +
            sourcevalues[vertexindexs[currentindex]];
    }
}

vector<float> lossdistgrid(const vector<float>& lossinputval,
                           int columncountv,
                           int rowcountvalu) {
    vector<float> distancemapv(columncountv * rowcountvalu, 999999.0f);

    cout << "    computing distance to loss pixels...\n";

    vector<unsigned char> losspixelmap(columncountv * rowcountvalu, 0);
    long long losspixelcnt = 0;

    for(int pixelindexer = 0;
         pixelindexer < columncountv * rowcountvalu;
         ++pixelindexer) {
        if(lossinputval[pixelindexer] >= 19.0f &&
            lossinputval[pixelindexer] <= 23.0f) {
            losspixelmap[pixelindexer] = 1;
            ++losspixelcnt;
        }
    }

    cout << "    found " << losspixelcnt << " loss pixels\n";

    if(losspixelcnt == 0) {
        return distancemapv;
    }

    const float infinitecost = 1e20f;

    vector<float> tempdistgrid(columncountv * rowcountvalu, infinitecost);

    const int maxdimension = max(columncountv, rowcountvalu);
    vector<int>   vertexindexv(maxdimension);
    vector<float> breakpgridvx(maxdimension + 1);
    vector<float> sourcegridvx(maxdimension);
    vector<float> distgridvalu(maxdimension);

    int barwidthvalu = 50;

    for(int columnindexv = 0; columnindexv < columncountv; ++columnindexv) {
        for(int rowindexvalu = 0;
             rowindexvalu < rowcountvalu;
             ++rowindexvalu) {
            int pixelindexer =
                rowindexvalu * columncountv + columnindexv;

            sourcegridvx[rowindexvalu] =
                losspixelmap[pixelindexer] ? 0.0f : infinitecost;
        }

        disttransone(sourcegridvx.data(),
                     rowcountvalu,
                     distgridvalu.data(),
                     vertexindexv.data(),
                     breakpgridvx.data());

        for(int rowindexvalu = 0;
             rowindexvalu < rowcountvalu;
             ++rowindexvalu) {
            tempdistgrid[rowindexvalu * columncountv + columnindexv] =
                distgridvalu[rowindexvalu];
        }

        float progressvalu =
            (float)(columnindexv + 1) / (float)columncountv;
        int   barfillvalu =
            (int)(barwidthvalu * progressvalu);

        cout << "    [";
        for(int barindexval = 0;
             barindexval < barwidthvalu;
             ++barindexval) {
            if(barindexval < barfillvalu) {
                cout << "=";
            } else if(barindexval == barfillvalu) {
                cout << ">";
            } else {
                cout << " ";
            }
        }
        cout << "] "
             << setw(3)
             << (int)(progressvalu * 100.0f)
             << "%\r";
        cout.flush();
    }

    cout << "\n";

    for(int rowindexvalu = 0;
         rowindexvalu < rowcountvalu;
         ++rowindexvalu) {
        int rowoffsetval = rowindexvalu * columncountv;

        for(int columnindexv = 0; columnindexv < columncountv; ++columnindexv) {
            sourcegridvx[columnindexv] =
                tempdistgrid[rowoffsetval + columnindexv];
        }

        disttransone(sourcegridvx.data(),
                     columncountv,
                     distgridvalu.data(),
                     vertexindexv.data(),
                     breakpgridvx.data());

        for(int columnindexv = 0; columnindexv < columncountv; ++columnindexv) {
            float distancepixel = sqrt(distgridvalu[columnindexv]);
            distancemapv[rowoffsetval + columnindexv] =
                distancepixel * 30.0f;
        }
    }

    return distancemapv;
}

vector<float> edgeforestid(const vector<float>& treeinputval,
                           const vector<float>& lossinputval,
                           int columncountv,
                           int rowcountvalu) {
    vector<float> edgemapvalue(columncountv * rowcountvalu, 0.0f);

    cout << "    detecting forest edges...\n";

    for(int rowindexvalu = 1;
         rowindexvalu < rowcountvalu - 1;
         ++rowindexvalu) {
        for(int columnindexv = 1;
             columnindexv < columncountv - 1;
             ++columnindexv) {
            int pixelindexer =
                rowindexvalu * columncountv + columnindexv;

            if(treeinputval[pixelindexer] > 50.0f &&
                lossinputval[pixelindexer] == 0.0f) {
                int forestneighb = 0;

                for(int neighrowoffs = -1;
                     neighrowoffs <= 1;
                     ++neighrowoffs) {
                    for(int neighcoloffs = -1;
                         neighcoloffs <= 1;
                         ++neighcoloffs) {
                        if(neighrowoffs == 0 &&
                            neighcoloffs == 0) {
                            continue;
                        }

                        int neighindexer =
                            (rowindexvalu + neighrowoffs) *
                                columncountv +
                            (columnindexv + neighcoloffs);

                        if(treeinputval[neighindexer] > 50.0f &&
                            lossinputval[neighindexer] == 0.0f) {
                            ++forestneighb;
                        }
                    }
                }

                edgemapvalue[pixelindexer] =
                    (8.0f - forestneighb) / 8.0f;
            }
        }
    }

    return edgemapvalue;
}

int main() {
    cout << "================================================================\n";
    cout << "ani spatial vulnerability index calculator\n";
    cout << "innovation: predictive conservation prioritization\n";
    cout << "================================================================\n";

    int columncountv;
    int rowcountvalu;

    // read config from processed folder
    ifstream configstream("data/processed/config.txt");
    if(!configstream) {
        cerr << "error: run prep first\n";
        return 1;
    }
    configstream >> columncountv >> rowcountvalu;

    cout << "\n[1/5] loading "
         << columncountv << " x "
         << rowcountvalu << " grid...\n";

    auto loadbinarygf = [&](const string& gridnamestrx) {
        vector<float> gridvaluemap(columncountv * rowcountvalu);

        ifstream gridfileifss("data/processed/" + gridnamestrx + ".bin",
                              ios::binary);
        if(!gridfileifss) {
            cerr << "missing: " << gridnamestrx << ".bin\n";
            exit(1);
        }

        gridfileifss.read(
            reinterpret_cast<char*>(gridvaluemap.data()),
            gridvaluemap.size() * sizeof(float)
        );

        return gridvaluemap;
    };

    vector<float> lossgridvalu = loadbinarygf("loss");
    vector<float> carbongridxv = loadbinarygf("bio");
    vector<float> slopegridval = loadbinarygf("slope");
    vector<float> treegridvalu = loadbinarygf("tree");

    cout << "  ok datasets loaded\n";

    cout << "\n[2/5] computing vulnerability components...\n";

    cout << "  [a] normalizing carbon density...\n";
    vector<float> carbonindexg =
        normalizgrid(carbongridxv, 0.0f, 400.0f);

    cout << "  [b] normalizing slope...\n";
    vector<float> slopeindexgr =
        normalizgrid(slopegridval, 0.0f, 45.0f);

    cout << "  [c] computing proximity to loss...\n";
    vector<float> lossproxgrid =
        lossdistgrid(lossgridvalu, columncountv, rowcountvalu);

    for(float& distancevalu : lossproxgrid) {
        distancevalu =
            max(0.0f, 1.0f - (distancevalu / 1000.0f));
    }

    cout << "  [d] detecting forest edges...\n";
    vector<float> edgeindexgrd =
        edgeforestid(treegridvalu, lossgridvalu,
                     columncountv, rowcountvalu);

    cout << "\n[3/5] computing composite vulnerability index...\n";

    vector<float> vulnindexgrd(columncountv * rowcountvalu, 0.0f);
    long long forestintact = 0;

    for(size_t pixelindexer = 0;
         pixelindexer < vulnindexgrd.size();
         ++pixelindexer) {
        if(treegridvalu[pixelindexer] > 50.0f &&
            lossgridvalu[pixelindexer] == 0.0f) {
            ++forestintact;

            vulnindexgrd[pixelindexer] =
                weightcarbon * carbonindexg[pixelindexer] +
                weightsloper * slopeindexgr[pixelindexer] +
                weightproxim * lossproxgrid[pixelindexer] +
                weightedgext * edgeindexgrd[pixelindexer];
        }
    }

    cout << "  ok analyzed " << forestintact
         << " intact forest pixels\n";

    cout << "\n[4/5] identifying conservation priorities...\n";

    vector<float> scorevaluegv;
    scorevaluegv.reserve((size_t)forestintact);

    for(float scorevaluenw : vulnindexgrd) {
        if(scorevaluenw > 0.0f) {
            scorevaluegv.push_back(scorevaluenw);
        }
    }

    auto compthresval = [&](double fractionvalu) -> float {
        size_t threshindexv =
            static_cast<size_t>(
                scorevaluegv.size() * fractionvalu
            );

        if(threshindexv >= scorevaluegv.size()) {
            threshindexv = scorevaluegv.size() - 1;
        }

        nth_element(scorevaluegv.begin(),
                    scorevaluegv.begin() + threshindexv,
                    scorevaluegv.end(),
                    greater<float>());

        return scorevaluegv[threshindexv];
    };

    float critthresold = compthresval(0.10);
    float highthresold = compthresval(0.25);
    float modrthresold = compthresval(0.50);

    long long critpixcount = 0;
    long long highpixcount = 0;
    long long modrpixcount = 0;

    double critcarbtotl = 0.0;
    double highcarbtotl = 0.0;
    double modrcarbtotl = 0.0;

    for(size_t pixelindexer = 0;
         pixelindexer < vulnindexgrd.size();
         ++pixelindexer) {
        float scorevaluenw = vulnindexgrd[pixelindexer];

        if(scorevaluenw >= critthresold) {
            ++critpixcount;
            critcarbtotl +=
                carbongridxv[pixelindexer] *
                pixelareacon * 0.47;
        } else if(scorevaluenw >= highthresold) {
            ++highpixcount;
            highcarbtotl +=
                carbongridxv[pixelindexer] *
                pixelareacon * 0.47;
        } else if(scorevaluenw >= modrthresold) {
            ++modrpixcount;
            modrcarbtotl +=
                carbongridxv[pixelindexer] *
                pixelareacon * 0.47;
        }
    }

    cout << "\n================================================================\n";
    cout << "conservation priority analysis\n";
    cout << "================================================================\n";

    cout << "\nhigh vulnerability forests (conservation priorities):\n";
    cout << "  +--------------+------------+--------------+--------------+\n";
    cout << "  | priority     |  area (ha) | carbon (t c) |  % of total  |\n";
    cout << "  +--------------+------------+--------------+--------------+\n";

    double foresthectar = forestintact * pixelareacon;

    cout << fixed << setprecision(1);

    cout << "  | critical     | "
         << setw(10) << (critpixcount * pixelareacon)
         << " | " << setw(12) << critcarbtotl
         << " | " << setw(11)
         << (100.0 * critpixcount / forestintact)
         << "% |\n";

    cout << "  | high         | "
         << setw(10) << (highpixcount * pixelareacon)
         << " | " << setw(12) << highcarbtotl
         << " | " << setw(11)
         << (100.0 * highpixcount / forestintact)
         << "% |\n";

    cout << "  | moderate     | "
         << setw(10) << (modrpixcount * pixelareacon)
         << " | " << setw(12) << modrcarbtotl
         << " | " << setw(11)
         << (100.0 * modrpixcount / forestintact)
         << "% |\n";

    cout << "  +--------------+------------+--------------+--------------+\n";

    cout << "\ninterpretation:\n";
    cout << "  - critical areas should receive immediate protection\n";
    cout << "  - high areas warrant monitoring and buffer zone creation\n";
    cout << "  - moderate areas suitable forsustainable use planning\n";

    cout << "\n[5/5] exporting vulnerability map...\n";

    ofstream binaryoutfil("data/processed/vulnerability-index.bin",
                          ios::binary);
    binaryoutfil.write(
        reinterpret_cast<char*>(vulnindexgrd.data()),
        vulnindexgrd.size() * sizeof(float)
    );
    binaryoutfil.close();

    cout << "  ok created: data/processed/vulnerability-index.bin\n";
    cout << "  (load in qgis to visualize priority zones)\n";

    ofstream csvoutputfil("vulnerability-summary.csv");
    csvoutputfil
        << "priorityclass,areahectares,carbontonnes,percentintact\n";
    csvoutputfil
        << "critical,"
        << (critpixcount * pixelareacon) << ","
        << critcarbtotl << ","
        << (100.0 * critpixcount / forestintact)
        << "\n";
    csvoutputfil
        << "high,"
        << (highpixcount * pixelareacon) << ","
        << highcarbtotl << ","
        << (100.0 * highpixcount / forestintact)
        << "\n";
    csvoutputfil
        << "moderate,"
        << (modrpixcount * pixelareacon) << ","
        << modrcarbtotl << ","
        << (100.0 * modrpixcount / forestintact)
        << "\n";
    csvoutputfil.close();

    cout << "  ok created: vulnerability-summary.csv\n";

    cout << "\n================================================================\n";
    cout << "analysis complete - conservation priorities identified\n";
    cout << "================================================================\n";
    return 0;
}