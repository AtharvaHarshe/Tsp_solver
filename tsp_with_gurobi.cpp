#include "gurobi_c++.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

// ==========================================================
// Compute Euclidean distance
// ==========================================================
double dist(double x1, double y1, double x2, double y2) {
    double dx = x1 - x2;
    double dy = y1 - y2;
    return sqrt(dx*dx + dy*dy);
}

// ==========================================================
// Find the shortest subtour in a solution
// EXACT SAME LOGIC AS OFFICIAL GUROBI TSP EXAMPLE
// ==========================================================
vector<int> find_subtour(const vector<vector<double>>& sol) {
    int n = sol.size();
    vector<bool> seen(n, false);
    vector<int> tour, best_tour;

    for (int start = 0; start < n; start++) {
        if (!seen[start]) {
            tour.clear();
            int i = start;

            while (!seen[i]) {
                seen[i] = true;
                tour.push_back(i);

                // follow the single outgoing edge
                for (int j = 0; j < n; j++) {
                    if (sol[i][j] > 0.5 && i != j) {
                        i = j;
                        break;
                    }
                }
            }

            if (best_tour.empty() || tour.size() < best_tour.size())
                best_tour = tour;
        }
    }
    return best_tour;
}

// ==========================================================
// Callback for lazy subtour elimination
// ==========================================================
class SubtourCallback : public GRBCallback {
public:
    int n;
    vector<vector<GRBVar>>& x;

    SubtourCallback(int n, vector<vector<GRBVar>>& x) :
        n(n), x(x) {}

protected:
    void callback() {
        try {
            if (where == GRB_CB_MIPSOL) {

                // Extract the current integer solution
                vector<vector<double>> sol(n, vector<double>(n));
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        sol[i][j] = getSolution(x[i][j]);

                // Find the shortest subtour
                vector<int> subtour = find_subtour(sol);

                if (subtour.size() < n) {
                    GRBLinExpr expr = 0;
                    for (int i = 0; i < subtour.size(); i++) {
                        for (int j = i + 1; j < subtour.size(); j++) {
                            int a = subtour[i];
                            int b = subtour[j];
                            expr += x[a][b] + x[b][a];
                        }
                    }
                    addLazy(expr <= (int)subtour.size() - 1);
                }
            }
        } catch (...) {}
    }
};

// ==========================================================
// MAIN PROGRAM
// ==========================================================
int main() {
    try {
        string filename = "TSP_1000_euclidianDistance.txt";
        ifstream fin(filename);
        if (!fin) {
            cout << "Could not open file.\n";
            return 1;
        }

        int n;
        fin >> n;

        vector<double> X(n), Y(n);
        string temp1, temp2;
        int i, j;
        double d;

        // Skip header line "Node1 Node2 Distance"
        fin >> temp1 >> temp2 >> temp1;

        // Read coordinates implicitly: we don't actually have coordinates,
        // we will create a complete graph by storing edges
        vector<vector<double>> distmat(n, vector<double>(n, 0));

        while (fin >> i >> j >> d) {
            i--; j--;
            distmat[i][j] = d;
            distmat[j][i] = d;
        }

        cout << "Loaded TSP with " << n << " nodes.\n";

        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "tsp_gurobi.log");
        env.start();

        GRBModel model = GRBModel(env);
        model.set(GRB_IntParam_LazyConstraints, 1);

        // Variables: x[i][j] = 1 if edge i→j used
        vector<vector<GRBVar>> x(n, vector<GRBVar>(n));

        for (int a = 0; a < n; a++) {
            for (int b = 0; b < n; b++) {
                if (a == b)
                    x[a][b] = model.addVar(0, 0, 0, GRB_BINARY);
                else
                    x[a][b] = model.addVar(0, 1, distmat[a][b], GRB_BINARY);
            }
        }

        // Out-degree = 1
        for (int a = 0; a < n; a++) {
            GRBLinExpr row = 0;
            for (int b = 0; b < n; b++) row += x[a][b];
            model.addConstr(row == 1);
        }

        // In-degree = 1
        for (int b = 0; b < n; b++) {
            GRBLinExpr col = 0;
            for (int a = 0; a < n; a++) col += x[a][b];
            model.addConstr(col == 1);
        }

        // Attach callback
        SubtourCallback cb(n, x);
        model.setCallback(&cb);

        model.set(GRB_DoubleParam_TimeLimit, 60.0);   // new in Gurobi 12
model.set(GRB_IntParam_MIPFocus, 1);
model.set(GRB_DoubleParam_Heuristics, 0.2);


        cout << "Solving...\n";
        model.optimize();

        if (model.get(GRB_IntAttr_SolCount) == 0) {
            cout << "No tour found within time.\n";
            return 0;
        }

        // Extract tour from solution
        vector<int> tour, seen(n, 0);
        int curr = 0;
        tour.push_back(curr);

        for (int k = 0; k < n - 1; k++) {
            for (int j = 0; j < n; j++) {
                if (x[curr][j].get(GRB_DoubleAttr_X) > 0.5) {
                    curr = j;
                    tour.push_back(curr);
                    break;
                }
            }
        }

        cout << "Best tour length: " 
             << model.get(GRB_DoubleAttr_ObjVal) << "\nTour:\n";

        for (int v : tour) cout << (v + 1) << ",";
        cout << "\n";

    } catch (GRBException& e) {
        cout << "Gurobi error: " << e.getMessage() << endl;
    }
}
