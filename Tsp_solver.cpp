#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <random>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <bits/stdc++.h>

using namespace std;
using clock_tp = chrono::high_resolution_clock;
using Clock = chrono::high_resolution_clock;

//Ecludian solver with its functions---------------






// ---------------- compute tour length ----------------
double tourLen(const vector<int>& T, const vector<vector<double>>& dist) {
    double L = 0.0;
    int n = (int)T.size();
    for (int i = 0; i < n - 1; i++)
        L += dist[T[i]][T[i + 1]];
    L += dist[T[n - 1]][T[0]];
    return L;
}

// ---------------- Best-of-N Nearest Neighbor ----------------
vector<int> bestNN(const vector<vector<double>>& dist,
                   mt19937_64& rng,
                   int trials = 40)
{
    int n = (int)dist.size();
    vector<int> bestTour;
    double bestCost = 1e18;
    uniform_int_distribution<int> pickStart(0, n - 1);

    for (int t = 0; t < trials; t++) {
        int start = pickStart(rng);
        vector<bool> used(n, false);
        vector<int> T;
        T.reserve(n);

        int cur = start;
        used[cur] = true;
        T.push_back(cur);

        for (int k = 1; k < n; k++) {
            double best = 1e18;
            int nxt = -1;
            for (int j = 0; j < n; j++) {
                if (!used[j] && dist[cur][j] < best) {
                    best = dist[cur][j];
                    nxt = j;
                }
            }
            cur = nxt;
            used[cur] = true;
            T.push_back(cur);
        }

        double L = tourLen(T, dist);
        if (L < bestCost) {
            bestCost = L;
            bestTour = T;
        }
    }

    return bestTour;
}

// ------------- Build candidate lists (K nearest neighbors) -------------
vector<vector<int>> buildCandidate(const vector<vector<double>>& dist, int K = 50) {
    int n = (int)dist.size();
    vector<vector<int>> cand(n);

    for (int i = 0; i < n; i++) {
        vector<pair<double, int>> tmp;
        tmp.reserve(n - 1);
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            tmp.push_back({ dist[i][j], j });
        }
        sort(tmp.begin(), tmp.end());
        int useK = min(K, (int)tmp.size());
        cand[i].reserve(useK);
        for (int k = 0; k < useK; k++)
            cand[i].push_back(tmp[k].second);
    }
    return cand;
}

// ---------------- FAST 2-opt (Candidate list based) ----------------
bool twoOptFast(vector<int>& T,
                const vector<vector<double>>& dist,
                const vector<vector<int>>& cand)
{
    int n = (int)T.size();
    bool improved = false;

    // position of each node in tour
    vector<int> pos(n);
    for (int i = 0; i < n; i++) pos[T[i]] = i;

    for (int idx = 0; idx < n; idx++) {
        int a = T[idx];
        int idxNext = (idx + 1 == n) ? 0 : (idx + 1);
        int b = T[idxNext];

        for (int nb : cand[a]) {
            int j = pos[nb];
            int jNext = (j + 1 == n) ? 0 : (j + 1);
            int d = T[jNext];

            if (idx == j || idxNext == j) continue;

            double oldCost = dist[a][b] + dist[nb][d];
            double newCost = dist[a][nb] + dist[b][d];

            if (newCost + 1e-12 < oldCost) {
                // Perform 2-opt reversal between (idxNext .. j)
                if (idxNext < j) {
                    reverse(T.begin() + idxNext, T.begin() + j + 1);
                } else {
                    reverse(T.begin() + jNext, T.begin() + idx + 1);
                }

                // Update positions
                for (int k = 0; k < n; k++)
                    pos[T[k]] = k;

                improved = true;
            }
        }
    }

    return improved;
}

// ------------- Or-opt limited (1- and 2-node relocations) -------------
bool orOptLimited(vector<int>& T,
                  const vector<vector<double>>& dist,
                  int MAX_MOVES = 2000)
{
    int n = (int)T.size();
    bool improved = false;
    int moves = 0;

    for (int seg = 1; seg <= 2; seg++) {
        for (int i = 0; i < n && moves < MAX_MOVES; i++) {
            int segEnd = i + seg - 1;
            if (segEnd >= n) break;

            for (int j = 0; j < n && moves < MAX_MOVES; j++) {
                if (j >= i && j <= segEnd) continue;

                int a = (i - 1 + n) % n;
                int b = (segEnd + 1) % n;
                int A = j;
                int C = (j + 1) % n;

                double oldCost = dist[T[a]][T[i]] +
                                 dist[T[segEnd]][T[b]] +
                                 dist[T[A]][T[C]];

                double newCost = dist[T[a]][T[b]] +
                                 dist[T[A]][T[i]] +
                                 dist[T[segEnd]][T[C]];

                if (newCost + 1e-12 < oldCost) {
                    vector<int> newT;
                    newT.reserve(n);

                    for (int k = 0; k < n; k++) {
                        if (k == (A + 1) % n) {
                            for (int t = i; t <= segEnd; t++)
                                newT.push_back(T[t]);
                        }
                        if (k < i || k > segEnd)
                            newT.push_back(T[k]);
                    }

                    T.swap(newT);
                    improved = true;
                    moves++;
                }
            }
        }
    }

    return improved;
}

// ---------------- Strong Double Bridge Perturbation ----------------
void doubleBridgeStrong(vector<int>& T, mt19937_64& rng) {
    int n = (int)T.size();
    if (n < 8) return;

    uniform_int_distribution<int> dist(0, n - 1);
    int a = dist(rng), b = dist(rng), c = dist(rng), d = dist(rng);
    vector<int> idx = { a,b,c,d };
    sort(idx.begin(), idx.end());
    a = idx[0]; b = idx[1]; c = idx[2]; d = idx[3];

    vector<int> p1(T.begin(), T.begin() + a);
    vector<int> p2(T.begin() + a, T.begin() + b);
    vector<int> p3(T.begin() + b, T.begin() + c);
    vector<int> p4(T.begin() + c, T.begin() + d);
    vector<int> p5(T.begin() + d, T.end());

    vector<int> nt;
    nt.reserve(n);

    // reorder segments to make a big "kick"
    nt.insert(nt.end(), p1.begin(), p1.end());
    nt.insert(nt.end(), p3.begin(), p3.end());
    nt.insert(nt.end(), p2.begin(), p2.end());
    nt.insert(nt.end(), p4.begin(), p4.end());
    nt.insert(nt.end(), p5.begin(), p5.end());

    T.swap(nt);
}

//main euclidian solver
vector<int> runEuclideanSolver(string fileName){
    auto global_start = clock_tp::now();  // measure full program time
    cout << "Disclamer: Works only for Eculidian distsnces (random distances may give incorrect answer)\n";
    

    // RNG
    mt19937_64 rng(
        (uint64_t)chrono::steady_clock::now().time_since_epoch().count()
    );

    // ---------------- read input file name ----------------
    

    string fname = fileName;
    cout << "Reading file: " << fname << "\n";

    ifstream f(fname);
if (!f) {
    cerr << "ERROR: cannot open file " << fname << "\n";
    return {};   // empty tour
}

int N;
if (!(f >> N)) {
    cerr << "ERROR: failed to read N from " << fname << "\n";
    return {};
}

if (N <= 0 || N > 20000) {
    cerr << "ERROR: N out of range (" << N << ") in " << fname << "\n";
    return {};
}

string h1, h2, h3;
if (!(f >> h1 >> h2 >> h3)) {
    cerr << "ERROR: failed to read header in " << fname << "\n";
    return {};
}

   

    
    

    vector<vector<double>> dist(N, vector<double>(N, 0.0));
    int u, v;
    double w;

    while (f >> u >> v >> w) {
        u--; v--;
        dist[u][v] = w;
        dist[v][u] = w;
    }
    f.close();
    
    // -------- Build candidate lists (once) --------
    auto candList = buildCandidate(dist, 50);

    // -------- Initial tour: best-of-40 NN + local search --------
    vector<int> tour = bestNN(dist, rng, 40);
    twoOptFast(tour, dist, candList);
    orOptLimited(tour, dist);
    twoOptFast(tour, dist, candList);

    double bestL = tourLen(tour, dist);
    vector<int> best = tour;

    cout << "Initial tour length: " << bestL << "\n";

    // -------- ILS Loop (hard time limit ~55s) --------
    const double LIMIT = 55.0;  // seconds for search part
    auto search_start = clock_tp::now();
    int iter = 0;

    while (true) {
        double t = chrono::duration<double>(clock_tp::now() - search_start).count();
        if (t > LIMIT) break;

        vector<int> T;

        // 70% of the time: perturb current best
        // 30% of the time: new NN start (to diversify more)    
        uniform_real_distribution<double> prob(0.0, 1.0);
        if (prob(rng) < 0.7) {
            T = best;
            doubleBridgeStrong(T, rng);
        } else {
            T = bestNN(dist, rng, 5);  // smaller trials for speed
        }

        // Local search
        twoOptFast(T, dist, candList);
        orOptLimited(T, dist);
        twoOptFast(T, dist, candList);

        double L = tourLen(T, dist);

        if (L < bestL) {
            bestL = L;
            best = T;
            cout << "Improved: " << L << "  at t=" << t << " sec\n";
        }

        iter++;
    }
    
    double tot = chrono::duration<double>(clock_tp::now() - search_start).count();
    cout << "time taken = " << tot << "\n" ;
    cout << "cost = " << bestL << "\n";
   
    return best;
    
}

//----------------------

//Random Solver---------------

// --------------------- Tour length ---------------------
double tourLength(const vector<int>& tour,
                  const vector<vector<double>>& dist) {
    double L = 0.0;
    int n = (int)tour.size();
    for (int i = 0; i < n - 1; ++i)
        L += dist[tour[i]][tour[i+1]];
    L += dist[tour[n-1]][tour[0]];
    return L;
}

// --------------------- Build full distance matrix ---------------------
// Input format:
// N
// Node1 Node2 Distance
// Node1 Node2 Distance
// ...
bool loadInstance(const string& fname,
                  vector<vector<double>>& dist) {
    ifstream fin(fname);
    if (!fin) {
        cerr << "Error: cannot open file " << fname << "\n";
        return false;
    }

    int N;
    fin >> N;
    if (!fin) {
        cerr << "Error: failed to read N\n";
        return false;
    }

    // Read header line (e.g. "Node1 Node2 Distance"), but don't rely on it.
    string h1, h2, h3;
    fin >> h1 >> h2 >> h3;

    dist.assign(N, vector<double>(N, 0.0));

    int u, v;
    double w;
    while (fin >> u >> v >> w) {
        // convert 1-based to 0-based
        --u; --v;
        if (u < 0 || v < 0 || u >= N || v >= N) continue;
        dist[u][v] = w;
        dist[v][u] = w;
    }
    return true;
}

// --------------------- Candidate lists (by cheap edges) ---------------------
// For random weights, we build for each node a list of its K cheapest neighbors.
vector<vector<int>> buildCandidates(const vector<vector<double>>& dist,
                                    int K = 40) {
    int N = (int)dist.size();
    vector<vector<int>> cand(N);

    for (int i = 0; i < N; ++i) {
        vector<pair<double,int>> tmp;
        tmp.reserve(N - 1);
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            tmp.push_back({dist[i][j], j});
        }
        sort(tmp.begin(), tmp.end(),
             [](const auto& a, const auto& b){ return a.first < b.first; });
        int useK = min(K, (int)tmp.size());
        cand[i].reserve(useK);
        for (int k = 0; k < useK; ++k)
            cand[i].push_back(tmp[k].second);
    }
    return cand;
}

// --------------------- 2-opt local search (candidate based) ---------------------
bool twoOptDescent(vector<int>& tour,
                   const vector<vector<double>>& dist,
                   const vector<vector<int>>& cand,
                   const Clock::time_point& deadline)
{
    int N = (int)tour.size();
    vector<int> pos(N);
    for (int i = 0; i < N; ++i) pos[tour[i]] = i;

    bool improved = false;

    while (true) {
        bool anyImproved = false;

        for (int i = 0; i < N; ++i) {
            // Time check
            if (Clock::now() >= deadline)
                return improved; // stop gently

            int a = tour[i];
            int iNext = (i + 1 == N) ? 0 : (i + 1);
            int b = tour[iNext];

            const auto& candList = cand[a];
            for (int nb : candList) {
                int j = pos[nb];
                int jNext = (j + 1 == N) ? 0 : (j + 1);
                int c = tour[j];
                int d = tour[jNext];

                // Avoid adjacent edges & trivial reversals
                if (i == j || iNext == j) continue;
                if (j == i || jNext == i) continue;

                double oldCost = dist[a][b] + dist[c][d];
                double newCost = dist[a][c] + dist[b][d];
                if (newCost + 1e-12 < oldCost) {
                    // perform 2-opt between (iNext .. j)
                    if (iNext < j) {
                        reverse(tour.begin() + iNext, tour.begin() + j + 1);
                    } else {
                        // wrap-around case
                        reverse(tour.begin() + jNext, tour.begin() + i + 1);
                    }
                    // update positions
                    for (int k = 0; k < N; ++k)
                        pos[tour[k]] = k;
                    anyImproved = true;
                    improved = true;
                    break; // restart scanning i
                }
            }

            if (anyImproved)
                break;
        }

        if (!anyImproved)
            break;
    }

    return improved;
}

// --------------------- Random 3-opt sampling ---------------------
// Try a limited number of random 3-opt moves to escape local minima.
// If an improving move is found, apply it and return true.
bool random3OptKick(vector<int>& tour,
                    const vector<vector<double>>& dist,
                    mt19937_64& rng,
                    int maxTries,
                    const Clock::time_point& deadline)
{
    int N = (int)tour.size();
    if (N < 6) return false;

    uniform_int_distribution<int> pick(0, N - 1);

    // position array for quick index lookups
    vector<int> pos(N);
    for (int i = 0; i < N; ++i) pos[tour[i]] = i;

    for (int attempt = 0; attempt < maxTries; ++attempt) {
        if (Clock::now() >= deadline) return false;

        int i = pick(rng);
        int j = pick(rng);
        int k = pick(rng);
        // ensure distinct and ordered indices
        vector<int> idx = {i,j,k};
        sort(idx.begin(), idx.end());
        i = idx[0]; j = idx[1]; k = idx[2];

        // need at least spacing
        if (i == j || j == k) continue;
        int i1 = i;
        int j1 = j;
        int k1 = k;

        int a = tour[i1];
        int b = tour[(i1 + 1) % N];
        int c = tour[j1];
        int d = tour[(j1 + 1) % N];
        int e = tour[k1];
        int f = tour[(k1 + 1) % N];

        double base = dist[a][b] + dist[c][d] + dist[e][f];

        // We will test a few standard 3-opt reconnections.
        // For simplicity, just implement 3 variants.

        // Option 1: reconnect (a-c), (b-e), (d-f)
        double opt1 = dist[a][c] + dist[b][e] + dist[d][f];

        // Option 2: (a-d), (e-b), (c-f)
        double opt2 = dist[a][d] + dist[e][b] + dist[c][f];

        // Option 3: (a-e), (d-b), (c-f)
        double opt3 = dist[a][e] + dist[d][b] + dist[c][f];

        double bestDelta = 0.0;
        int bestType = 0;

        if (opt1 + 1e-12 < base) {
            bestDelta = opt1 - base;
            bestType = 1;
        }
        if (opt2 + 1e-12 < base && opt2 - base < bestDelta) {
            bestDelta = opt2 - base;
            bestType = 2;
        }
        if (opt3 + 1e-12 < base && opt3 - base < bestDelta) {
            bestDelta = opt3 - base;
            bestType = 3;
        }

        if (bestType == 0) continue; // no improvement

        // Apply chosen reconnection by rearranging segments
        // segments: [i1+1 .. j1], [j1+1 .. k1], etc.
        vector<int> newTour;
        newTour.reserve(N);

        auto addSegment = [&](int from, int to, bool rev) {
            // inclusive indices, possibly wrap-around
            vector<int> tmp;
            if (from <= to) {
                for (int x = from; x <= to; ++x)
                    tmp.push_back(tour[x]);
            } else {
                for (int x = from; x < N; ++x) tmp.push_back(tour[x]);
                for (int x = 0; x <= to; ++x) tmp.push_back(tour[x]);
            }
            if (rev) reverse(tmp.begin(), tmp.end());
            newTour.insert(newTour.end(), tmp.begin(), tmp.end());
        };

        // We'll only handle non-wrap segments assuming 0 <= i1 < j1 < k1 < N
        // which is true by construction. That simplifies implementation.
        //
        // segments in order: [0..i1], [i1+1..j1], [j1+1..k1], [k1+1..N-1]
        // We'll handle each option by choosing how to permute/reverse these segments.

        const int A0 = 0, A1 = i1;
        const int B0 = i1+1, B1 = j1;
        const int C0 = j1+1, C1 = k1;
        const int D0 = k1+1, D1 = N-1;

        newTour.clear();
        if (bestType == 1) {
            // Option 1 pattern:
            // Keep [0..i1], then reverse [B], then [C], then [D]
            addSegment(A0, A1, false);
            addSegment(B0, B1, true);
            addSegment(C0, C1, false);
            addSegment(D0, D1, false);
        } else if (bestType == 2) {
            // Option 2:
            // [0..i1], [C rev], [B], [D]
            addSegment(A0, A1, false);
            addSegment(C0, C1, true);
            addSegment(B0, B1, false);
            addSegment(D0, D1, false);
        } else if (bestType == 3) {
            // Option 3:
            // [0..i1], [C], [B rev], [D]
            addSegment(A0, A1, false);
            addSegment(C0, C1, false);
            addSegment(B0, B1, true);
            addSegment(D0, D1, false);
        }

        if ((int)newTour.size() == N) {
            tour.swap(newTour);
            return true;
        }
    }

    return false;
}

// --------------------- Generate random permutation tour ---------------------
vector<int> randomTour(int N, mt19937_64& rng) {
    vector<int> t(N);
    iota(t.begin(), t.end(), 0);
    shuffle(t.begin(), t.end(), rng);
    return t;
}


vector<int> runRandomSolver(string fileName, int timeLimitSeconds){
    auto global_start = Clock::now();
    cout << "Disclamer: Works only for random distances (Eculidian distsnces may give incorrect answer)\n";
    // Hard total limit: 55 seconds
    const double TOTAL_LIMIT = 55.0;
    auto deadline =
        global_start + chrono::duration_cast<Clock::duration>(chrono::duration<double>(TOTAL_LIMIT));

    

    string fname = fileName;
    cout << "Reading file: " << fname << "\n";

   
   

    
        ifstream f(fname);
        if (!f) {
        cerr << "ERROR: cannot open file " << fname << "\n";
        return {};   // empty tour
    }

    int N;
    if (!(f >> N)) {
        cerr << "ERROR: failed to read N from " << fname << "\n";
        return {};
    }

    if (N <= 0 || N > 20000) {
        cerr << "ERROR: N out of range (" << N << ") in " << fname << "\n";
        return {};
    }

    string h1, h2, h3;
    if (!(f >> h1 >> h2 >> h3)) {
        cerr << "ERROR: failed to read header in " << fname << "\n";
        return {};
    }


    vector<vector<double>> dist(N, vector<double>(N, 0.0));
    int u, v;
    double w;

    while (f >> u >> v >> w) {
        u--; v--;
        dist[u][v] = w;
        dist[v][u] = w;
    }
    f.close();


    // RNG
    mt19937_64 rng(
        (uint64_t)chrono::steady_clock::now().time_since_epoch().count()
    );

    // Build candidate lists once
    auto cand = buildCandidates(dist, 40);

    // ---- Initial solution: simple nearest-neighbor from random start ----
    // (We do a few NN starts and keep the best one before local search)
    int NN_TRIALS = 10;
    vector<int> bestTour;
    double bestCost = 1e18;

    for (int t = 0; t < NN_TRIALS; ++t) {
        if (Clock::now() >= deadline) break;

        uniform_int_distribution<int> pickStart(0, N-1);
        int start = pickStart(rng);

        vector<bool> used(N,false);
        vector<int> T; T.reserve(N);
        int cur = start;
        used[cur] = true;
        T.push_back(cur);

        for (int k = 1; k < N; ++k) {
            double best = 1e18;
            int nxt = -1;
            for (int j = 0; j < N; ++j) {
                if (!used[j] && dist[cur][j] < best) {
                    best = dist[cur][j];
                    nxt = j;
                }
            }
            cur = nxt;
            used[cur] = true;
            T.push_back(cur);
        }

        double cost = tourLength(T, dist);
        if (cost < bestCost) {
            bestCost = cost;
            bestTour = T;
        }
    }

    // Local search on this best NN tour
    if (!bestTour.empty()) {
        twoOptDescent(bestTour, dist, cand, deadline);
        // Optional few 3-opt kicks and re-2-opt
        for (int k = 0; k < 3; ++k) {
            if (Clock::now() >= deadline) break;
            if (random3OptKick(bestTour, dist, rng, 200, deadline)) {
                twoOptDescent(bestTour, dist, cand, deadline);
            } else break;
        }
        bestCost = tourLength(bestTour, dist);
    } else {
        // Fallback: random tour
        bestTour = randomTour(N, rng);
        bestCost = tourLength(bestTour, dist);
    }

    cerr << "Initial best cost: " << bestCost << "\n";

    // ---- Multi-start loop ----
    int iterations = 0;
    while (Clock::now() < deadline) {
        vector<int> T = randomTour(N, rng);

        twoOptDescent(T, dist, cand, deadline);
        if (Clock::now() >= deadline) break;

        // a couple of random 3-opt kicks
        for (int k = 0; k < 2; ++k) {
            if (Clock::now() >= deadline) break;
            if (random3OptKick(T, dist, rng, 300, deadline)) {
                twoOptDescent(T, dist, cand, deadline);
            } else {
                break;
            }
        }

        double cost = tourLength(T, dist);
        if (cost < bestCost) {
            bestCost = cost;
            bestTour = T;
            cerr << "Improved: " << bestCost
                 << " at t = "
                 << chrono::duration<double>(Clock::now() - global_start).count()
                 << " sec\n";
        }

        ++iterations;
    }

    double tot = chrono::duration<double>(Clock::now() - global_start).count();
    cout << "time taken = " << tot << "\n" ;
    cout << "cost = " << bestCost << "\n";

    
    return bestTour;
}



//--------------------------





//main code

int main() {
    string euclidFile = "TSP_1000_euclidianDistance.txt";
    string randomFile = "TSP_1000_randomDistance.txt";
    string outputFile = "solution_924610301.txt";

    int timeLimit = 55;     // seconds each solver runs

    cout << "Solving Euclidean Instance...\n";
    vector<int> tourE = runEuclideanSolver(euclidFile);

    // --- Write Euclidean result to file (overwrite) ---
    ofstream out(outputFile, ios::out);
    if(!out.is_open()) { 
        cerr<<"ERROR: Cannot create output file\n"; 
        return 0;
    }

    for(int i=0;i<tourE.size();i++){
        out << tourE[i]+1 << (i+1<tourE.size() ? "," : "");
    }
    out << "," << tourE[0]+1;  // to complete the cycle
    out << "\n";   // first line complete
    out.close();
    cout << "Euclidean solution stored on line 1.\n";


    cout << "Solving Random Instance...\n";
    vector<int> tourR = runRandomSolver(randomFile, timeLimit);

    // --- Append Random result into 2nd line ---
    out.open(outputFile, ios::app);
    for(int i=0;i<tourR.size();i++){
        out << tourR[i]+1 << (i+1<tourR.size() ? "," : "");
    }
    out << "," << tourR[0]+1;  // to complete the cycle
    out << "\n";
    out.close();
    cout << "Random solution appended on line 2.\n";


    cout << "Task Finished → output stored in: " << outputFile << endl;
    return 0;
}
