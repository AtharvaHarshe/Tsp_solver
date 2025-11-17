#include <bits/stdc++.h>
using namespace std;
using clock_tp = chrono::high_resolution_clock;

// -------------------- compute tour length --------------------
double tourLen(const vector<int>& T, const vector<vector<double>>& dist) {
    double L=0; int n=T.size();
    for(int i=0;i<n-1;i++) L+=dist[T[i]][T[i+1]];
    return L + dist[T[n-1]][T[0]];
}

// -------------------- Nearest Neighbor --------------------
vector<int> nearestNeighbor(const vector<vector<double>>& dist) {
    int n=dist.size();
    vector<bool> used(n,false);
    vector<int> tour; tour.reserve(n);
    int cur=0;
    used[cur]=true; tour.push_back(cur);

    for(int k=1;k<n;k++){
        double best=1e18; int nxt=-1;
        for(int j=0;j<n;j++){
            if(!used[j] && dist[cur][j]<best){
                best=dist[cur][j]; nxt=j;
            }
        }
        cur=nxt; used[cur]=true;
        tour.push_back(cur);
    }
    return tour;
}

// -------------------- FAST 2-opt --------------------
bool twoOpt(vector<int>& T, const vector<vector<double>>& dist) {
    int n=T.size();
    bool improved=false;

    for(int i=1;i<n-2;i++){
        for(int j=i+1;j<n-1;j++){
            double old = dist[T[i-1]][T[i]] + dist[T[j]][T[j+1]];
            double nw  = dist[T[i-1]][T[j]] + dist[T[i]][T[j+1]];

            if(nw < old){
                reverse(T.begin()+i, T.begin()+j+1);
                improved=true;
            }
        }
    }
    return improved;
}

// -------------------- FAST Or-opt limited to K nearest tries --------------------
bool orOptLimited(vector<int>& T, const vector<vector<double>>& dist, int MAX_MOVES=2000) {
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

                int A = j;
                int C = (j+1)%n;

                double old = dist[T[a]][T[i]] + dist[T[segEnd]][T[b]]
                           + dist[T[A]][T[C]];

                double now = dist[T[a]][T[b]] + dist[T[A]][T[i]]
                           + dist[T[segEnd]][T[C]];

                if(now < old){
                    // build new tour fast
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

// -------------------- Double Bridge --------------------
void doubleBridge(vector<int>& T){
    int n=T.size();
    vector<int> idx(n);
    iota(idx.begin(),idx.end(),0);
    shuffle(idx.begin(), idx.end(), mt19937_64(random_device{}()));
    sort(idx.begin(), idx.begin()+4);

    int a=idx[0], b=idx[1], c=idx[2], d=idx[3];

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

// -------------------- MAIN --------------------
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string fname="TSP_1000_euclidianDistance.txt";
    ifstream f(fname);
    if(!f){ cerr<<"Cannot open file\n"; return 1; }

    int N; f>>N;
    string h1,h2,h3; f>>h1>>h2>>h3;

    vector<vector<double>> dist(N, vector<double>(N,0));

    int u,v; double w;
    while(f>>u>>v>>w){
        u--; v--;
        dist[u][v]=w;
        dist[v][u]=w;
    }
    f.close();

    // -------- Initial tour --------
    vector<int> tour = nearestNeighbor(dist);
    twoOpt(tour, dist);
    orOptLimited(tour, dist);
    twoOpt(tour, dist);

    double bestL = tourLen(tour, dist);
    vector<int> best = tour;

    cout<<"Initial LS length: "<<bestL<<"\n";

    // -------- ILS Loop --------
    const double LIMIT=55.0; // hard 55-second limit
    auto t0 = clock_tp::now();
    int iter=0;

    while(true){
        double t = chrono::duration<double>(clock_tp::now()-t0).count();
        if(t > LIMIT) break;

        vector<int> T = best;
        doubleBridge(T);

        twoOpt(T, dist);
        orOptLimited(T, dist);
        twoOpt(T, dist);

        double L = tourLen(T, dist);

        if(L < bestL){
            bestL = L;
            best = T;
            cout<<"Improved: "<<L<<"   time="<<t<<" sec\n";
        }
        iter++;
    }

    double TEND = chrono::duration<double>(clock_tp::now()-t0).count();

    cout<<"\nBest length: "<<bestL<<"\n";
    cout<<"Time: "<<TEND<<" sec\n";
    cout<<"Iterations: "<<iter<<"\n\nTour:\n";

    for(int i=0;i<N;i++) cout<<best[i]+1<<",";
    cout<<"\n";
}
