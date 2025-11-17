#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <random>
#include <chrono>
#include <cmath>
using namespace std;
using clock_tp = chrono::high_resolution_clock;

// ---------------- compute tour length ----------------
double tourLen(const vector<int>& T, const vector<vector<double>>& dist) {
    double L=0; 
    int n=T.size();
    for(int i=0;i<n-1;i++) L+=dist[T[i]][T[i+1]];
    return L + dist[T[n-1]][T[0]];
}

// ---------------- Best-of-40 Nearest Neighbor ----------------
vector<int> bestNN(const vector<vector<double>>& dist, int trials=40) {
    int n = dist.size();
    vector<int> bestTour;
    double bestCost = 1e18;

    for(int t=0; t<trials; t++){
        int start = rand() % n;
        vector<bool> used(n,false);
        vector<int> T; T.reserve(n);

        int cur = start;
        used[cur] = true;
        T.push_back(cur);

        for(int k=1;k<n;k++){
            double best = 1e18; 
            int nxt = -1;
            for(int j=0;j<n;j++){
                if(!used[j] && dist[cur][j] < best){
                    best = dist[cur][j]; 
                    nxt = j;
                }
            }
            cur = nxt;
            used[cur] = true;
            T.push_back(cur);
        }

        double L = tourLen(T,dist);
        if(L < bestCost){
            bestCost = L;
            bestTour = T;
        }
    }

    return bestTour;
}

// ------------- Build candidate lists (K nearest neighbors) -------------
vector<vector<int>> buildCandidate(const vector<vector<double>>& dist, int K=30) {
    int n = dist.size();
    vector<vector<int>> cand(n);

    for(int i=0;i<n;i++){
        vector<pair<double,int>> tmp;
        tmp.reserve(n-1);
        for(int j=0;j<n;j++){
            if(i == j) continue;
            tmp.push_back({dist[i][j], j});
        }
        sort(tmp.begin(), tmp.end());
        cand[i].reserve(K);
        for(int k=0;k<K;k++) cand[i].push_back(tmp[k].second);
    }
    return cand;
}

// ---------------- FAST 2-opt (Candidate list based) ----------------
bool twoOptFast(vector<int>& T, const vector<vector<double>>& dist,
                const vector<vector<int>>& cand)
{
    int n = T.size();
    bool improved = false;

    vector<int> pos(n);
    for(int i=0;i<n;i++) pos[T[i]] = i;

    for(int idx=0; idx<n; idx++){
        int a = T[idx];
        int a_next = T[(idx+1)%n];

        for(int b : cand[a]) {
            int j = pos[b];
            int b_next = T[(j+1)%n];

            double old = dist[a][a_next] + dist[b][b_next];
            double nw  = dist[a][b]      + dist[a_next][b_next];

            if(nw + 1e-12 < old){
                if(idx < j)
                    reverse(T.begin()+idx+1, T.begin()+j+1);
                else
                    reverse(T.begin()+j+1, T.begin()+idx+1);

                for(int k=0;k<n;k++) pos[T[k]] = k;
                improved = true;
            }
        }
    }
    return improved;
}

// ------------- Or-opt limited (1- and 2-node relocations) -------------
bool orOptLimited(vector<int>& T, const vector<vector<double>>& dist, int MAX_MOVES=1500) {
    int n=T.size();
    bool improved=false;
    int moves=0;

    for(int seg=1;seg<=2;seg++){
        for(int i=0;i<n && moves<MAX_MOVES;i++){
            int segEnd = i+seg-1;
            if(segEnd>=n) break;

            for(int j=0;j<n && moves<MAX_MOVES;j++){
                if(j>=i && j<=segEnd) continue;

                int a = (i-1+n)%n;
                int b = (segEnd+1)%n;

                int A=j;
                int C=(j+1)%n;

                double old = dist[T[a]][T[i]] + dist[T[segEnd]][T[b]]
                           + dist[T[A]][T[C]];

                double nw  = dist[T[a]][T[b]] + dist[T[A]][T[i]]
                           + dist[T[segEnd]][T[C]];

                if(nw + 1e-12 < old){
                    vector<int> newT; newT.reserve(n);

                    for(int k=0;k<n;k++){
                        if(k==(A+1)%n){
                            for(int t=i;t<=segEnd;t++) newT.push_back(T[t]);
                        }
                        if(k<i || k>segEnd) newT.push_back(T[k]);
                    }

                    T.swap(newT);
                    improved=true;
                    moves++;
                }
            }
        }
    }
    return improved;
}

// ---------------- Double Bridge Perturbation ----------------
void doubleBridge(vector<int>& T) {
    int n=T.size();
    vector<int> idx(n);
    iota(idx.begin(), idx.end(), 0);
    shuffle(idx.begin(), idx.end(), mt19937_64(random_device{}()));
    sort(idx.begin(), idx.begin()+4);

    int a = idx[0], b = idx[1], c = idx[2], d = idx[3];

    vector<int> p1(T.begin(), T.begin()+a);
    vector<int> p2(T.begin()+a, T.begin()+b);
    vector<int> p3(T.begin()+b, T.begin()+c);
    vector<int> p4(T.begin()+c, T.begin()+d);
    vector<int> p5(T.begin()+d, T.end());

    vector<int> nt;
    nt.reserve(n);
    nt.insert(nt.end(), p1.begin(), p1.end());
    nt.insert(nt.end(), p3.begin(), p3.end());
    nt.insert(nt.end(), p2.begin(), p2.end());
    nt.insert(nt.end(), p4.begin(), p4.end());
    nt.insert(nt.end(), p5.begin(), p5.end());

    T.swap(nt);
}

// =====================================================================
//                                MAIN
// =====================================================================
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string fname="TSP_1000_euclidianDistance.txt";
    ifstream f(fname);
    if(!f){ cerr<<"Cannot open file\n"; return 1; }

    int N;
    f >> N;
    string h1,h2,h3;
    f >> h1 >> h2 >> h3;

    vector<vector<double>> dist(N, vector<double>(N,0));
    int u,v; 
    double w;
    while(f >> u >> v >> w){
        u--; v--;
        dist[u][v] = w;
        dist[v][u] = w;
    }
    f.close();

    // -------- Build candidate lists (once) --------
    auto candList = buildCandidate(dist, 30);

    // -------- Best-of-40 NN initial tour --------
    vector<int> tour = bestNN(dist);
    twoOptFast(tour, dist, candList);
    orOptLimited(tour, dist);
    twoOptFast(tour, dist, candList);

    double bestL = tourLen(tour, dist);
    vector<int> best = tour;

    cout << "Initial tour length: " << bestL << "\n";

    // -------- ILS Loop (55 seconds) --------
    const double LIMIT = 60*1;
    auto t0 = clock_tp::now();
    int iter=0;

    while(true){
        double t = chrono::duration<double>(clock_tp::now() - t0).count();
        if(t > LIMIT) break;

        vector<int> T = best;   // intensification
        doubleBridge(T);        // diversification

        twoOptFast(T, dist, candList);
        orOptLimited(T, dist);
        twoOptFast(T, dist, candList);

        double L = tourLen(T, dist);

        if(L < bestL){
            bestL = L;
            best = T;
            cout << "Improved: " << L << "  at t=" << t << " sec\n";
        }

        iter++;
    }

    double TEND = chrono::duration<double>(clock_tp::now() - t0).count();
    cout << "\nBest length: " << bestL << "\n";
    cout << "Time: " << TEND << " sec\n";
    cout << "Iterations: " << iter << "\n\nTour:\n";

    for(int i=0;i<N;i++)
        cout << best[i]+1 << ",";
    cout << "\n";
}
