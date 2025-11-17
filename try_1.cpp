#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int main() {
    string filename = "TSP_1000_euclidianDistance.txt";   // your txt file name
    ifstream file(filename);

    if (!file.is_open()) {
        cout << "Error: Could not open file!" << endl;
        return 1;
    }

    // Read number of edges
    int N;
    file >> N;

    // Skip header (Node1 Node2 Distance)
    string h1, h2, h3;
    file >> h1 >> h2 >> h3;

    vector<int> node1(N), node2(N);
    vector<double> dist(N);

    for (int i = 0; i < N; i++) {
        file >> node1[i] >> node2[i] >> dist[i];
    }

    file.close();

    // Print some loaded values to check
    cout << "Loaded " << N << " edges." << endl;
    cout << "First 5 entries:\n";

    for (int i = 0; i < 50 && i < N; i++) {
        cout << node1[i] << "  " 
             << node2[i] << "  " 
             << dist[i]  << endl;
    }

    return 0;
}
