#include <bits/stdc++.h>
using namespace std;
using namespace std::chrono;

int main() {
    string filename = "TSP_1000_euclidianDistance.txt";
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file.\n";
        return 1;
    }

    int N;
    file >> N;  // Number of nodes = 1000

    string a, b, c;
    file >> a >> b >> c;  // "Node1 Node2 Distance"

    vector<vector<double>> dist(N, vector<double>(N, 0.0));

    int u, v;
    double w;

    // Read ALL 499,500 edges
    while (file >> u >> v >> w) {
        u--; v--;
        dist[u][v] = w;
        dist[v][u] = w;
    }

    // --------------------- NN TOUR ----------------------
    vector<int> tour;
    tour.reserve(N);
    vector<bool> visited(N,false);

    int current = 0;
    tour.push_back(current);
    visited[current] = true;

    for (int step = 1; step < N; step++) {
        double best = 1e18;
        int bestNode = -1;

        for (int j = 0; j < N; j++) {
            if (!visited[j] && dist[current][j] < best) {
                best = dist[current][j];
                bestNode = j;
            }
        }
        current = bestNode;
        visited[current] = true;
        tour.push_back(current);
    }

    // --------------------- 2-OPT ------------------------
    auto cost = [&](const vector<int>& T){
        double c=0;
        for(int i=0;i+1<T.size();i++)
            c += dist[T[i]][T[i+1]];
        c += dist[T.back()][T[0]];
        return c;
    };

    auto start = high_resolution_clock::now();

    bool improved = true;
    while (improved) {
        improved = false;
        for (int i = 1; i < N - 2; i++) {
            for (int j = i + 1; j < N - 1; j++) {
                double oldDist =
                    dist[tour[i - 1]][tour[i]] +
                    dist[tour[j]][tour[j + 1]];

                double newDist =
                    dist[tour[i - 1]][tour[j]] +
                    dist[tour[i]][tour[j + 1]];

                if (newDist < oldDist) {
                    reverse(tour.begin() + i, tour.begin() + j + 1);
                    improved = true;
                }
            }
        }
    }

    auto end = high_resolution_clock::now();
    double solve_time = duration<double>(end - start).count();

    double bestCost = cost(tour);

    // --------------------- OUTPUT ------------------------
    cout << "Best tour length: " << bestCost << endl;
    cout << "Solve time: " << solve_time << " seconds\n";
    cout << "Tour (1-based):\n";

    for (int x : tour) cout << x + 1 << ",";
    cout << "\n";

    return 0;
}
