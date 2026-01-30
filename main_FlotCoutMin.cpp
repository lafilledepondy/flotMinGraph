#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <fstream>
#include <stdexcept>
#include <assert.h>

using namespace std;

/***********************************************************************/
/************************ STRUCTURES DE DONNEES ************************/
/***********************************************************************/

struct Arc
{
    int toNode;
    int capacityResidu; // r_ij = MAXCAP (u_ij) - flow (x_ij)
    int MAXCAP;
    int MINCAP;
    int cost;         // c_ij for forward, -c_ij for reverse
    int reverseArcId; // Index of the corresponding reverse arc
};

struct Node
{
    int id;
    int supply;         // b_i for working supply (algorithm)
    int originalSupply; // Original supply from input file (for validation)
};

struct Graph
{
    int numNodes; // NODES
    int numArcs;  // DENSITY
    vector<Node> nodes;
    vector<vector<Arc>> adj; // Residual graph adjacency list
};

/***********************************************************************/
/************************** GRAPH UTILITIES ****************************/
/***********************************************************************/

Graph makeEmptyGraph(int n)
{
    Graph g;
    g.numNodes = n;
    g.numArcs = 0;
    g.nodes.resize(n);
    g.adj.resize(n);
    for (int i = 0; i < n; ++i)
    {
        g.nodes[i].id = i;
        g.nodes[i].supply = 0;
        g.nodes[i].originalSupply = 0;
    }
    return g;
}

bool isValidIndex(const Graph &g, int node_idx)
{
    return (node_idx >= 0 && node_idx < g.numNodes);
}

void addArc(Graph &g, int u, int v, int mincap, int maxcap, int cost)
{
    Arc fwd = {v, maxcap, maxcap, mincap, cost, -1};
    Arc rev = {u, 0, 0, 0, -cost, -1};

    g.adj[u].push_back(fwd);
    g.adj[v].push_back(rev);

    int fwdIdx = static_cast<int>(g.adj[u].size()) - 1;
    int revIdx = static_cast<int>(g.adj[v].size()) - 1;

    g.adj[u][fwdIdx].reverseArcId = revIdx;
    g.adj[v][revIdx].reverseArcId = fwdIdx;

    g.numArcs += 1; // Count only forward arcs
}

void printGraph(const Graph &g)
{
    cout << "Graph: " << g.numNodes << " nodes, " << g.numArcs << " arcs" << endl;
    for (int u = 0; u < g.numNodes; ++u)
    {
        cout << "Node " << u << " (orig_supply=" << g.nodes[u].originalSupply
             << "):" << endl;
        for (const Arc &arc : g.adj[u])
        {
            int flow = arc.MAXCAP - arc.capacityResidu;
            cout << "  -> " << arc.toNode
                 << " flow=" << flow
                 << "/" << arc.MAXCAP
                 << " cost=" << arc.cost
                 << " res_cap=" << arc.capacityResidu
                 << " rev_id=" << arc.reverseArcId << endl;
        }
    }
}

/***********************************************************************/
/****************************** READ FILE ******************************/
/***********************************************************************/

Graph readGraph(const string &filename)
{
    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Error: Cannot open the file" << endl;
        throw runtime_error("Cannot open input file");
    }

    // Skip comments
    while (true)
    {
        int c = file.peek();
        if (c == 'c' || c == '\n' || c == '\r' || c == EOF)
        {
            string dump;
            getline(file, dump);
        }
        else
        {
            break;
        }
    }

    // Read header: p min NODES DENSITY
    char p;
    string format;
    int numNodes, numArcs;
    file >> p >> format >> numNodes >> numArcs;
    if (!file.good() || p != 'p' || format != "min")
    {
        throw runtime_error("Invalid header format");
    }

    Graph g = makeEmptyGraph(numNodes);
    vector<int> supply(numNodes + 1, 0);

    // Read n and a lines
    char type;
    while (file >> type)
    {
        if (type == 'n')
        {
            int id, s;
            file >> id >> s;
            // Verify node id's bounds
            if (file.good() && id >= 1 && id <= numNodes)
            {
                supply[id] = s;
            }
        }
        else if (type == 'a')
        {
            int u, v, mincap, maxcap, cost;
            file >> u >> v >> mincap >> maxcap >> cost;
            // Verify arc endpoints bounds
            if (file.good() && u >= 1 && u <= numNodes && v >= 1 && v <= numNodes)
            {
                addArc(g, u - 1, v - 1, mincap, maxcap, cost);
            }
        }
        else
        {
            // Throw away unknown lines
            string dump;
            getline(file, dump);
        }
    }
    file.close();

    // Store supplies
    for (int i = 1; i <= numNodes; ++i)
    {
        g.nodes[i - 1].originalSupply = supply[i];
        g.nodes[i - 1].supply = supply[i];
    }

    return g;
}

/***********************************************************************/
/**************************** BELLMAN-FORD *****************************/
/***********************************************************************/

bool hasNegativeCycle(const Graph &g, const vector<int> &dist)
{
    const int n = g.numNodes;
    const int INF = numeric_limits<int>::max() / 4; // INF is a very large value.
                                                    // By /4 => making into a small value mainly to avoid overflow

    for (int u = 0; u < n; ++u)
    {
        if (dist[u] == INF)
        {
            continue; // la distance u is unreachable node
        }

        for (const Arc &arc : g.adj[u])
        {
            if (arc.capacityResidu <= 0)
            {
                continue; // only consider arcs with residual capacity
            }

            int v = arc.toNode;
            if (dist[u] + arc.cost < dist[v])
            {
                return true; // reachable negative cycle
            }
        }
    }
    return false;
}

bool updateDistances(const Graph &g, int u,
                     vector<int> &dist,
                     vector<int> &parentNode,
                     vector<int> &parentArc)
{
    bool updated = false;
    const int INF = numeric_limits<int>::max() / 4; // INF is a very large value.
                                                    // By /4 => making into a small value mainly to avoid overflow

    if (dist[u] == INF)
    {
        return false; // u is unreachable node
    }

    for (int idx = 0; idx < static_cast<int>(g.adj[u].size()); ++idx)
    {
        const Arc &arc = g.adj[u][idx];
        if (arc.capacityResidu <= 0)
        {
            continue; // only consider arcs with residual capacity
        }

        int v = arc.toNode;
        int ndist = dist[u] + arc.cost;

        if (ndist < dist[v]) // updating distances -> current distance is smaller than previously known
        {
            dist[v] = ndist;
            parentNode[v] = u;  // store predecessor node
            parentArc[v] = idx; // store predecessor arc
            updated = true;
        }
    }
    return updated;
}

bool bellmanFord(const Graph &g, int src,
                 vector<int> &dist,
                 vector<int> &parentNode,
                 vector<int> &parentArc)
{
    const int n = g.numNodes;
    const int INF = numeric_limits<int>::max() / 4;

    // Initialization
    dist.assign(n, INF);
    parentNode.assign(n, -1);
    parentArc.assign(n, -1);
    dist[src] = 0;

    // Relax edges n-1 times
    for (int k = 0; k < n - 1; ++k)
    {
        bool anyUpdated = false;
        for (int u = 0; u < n; ++u)
        {
            if (updateDistances(g, u, dist, parentNode, parentArc))
                anyUpdated = true;
        }
        if (!anyUpdated)
            break;
    }

    // Negative cycle detection
    if (hasNegativeCycle(g, dist))
        return false;

    return true;
}

/***********************************************************************/
/**************************** FLOT COUT MIN *****************************/
/***********************************************************************/

bool trouverCheminMin(const Graph &g, int s,
                      vector<int> &dist,
                      vector<int> &parentNode,
                      vector<int> &parentArc,
                      int &t)
{
    // si Bellman-Ford detecte un cycle negatif alors retourner false car pas de chemin de cout minimum
    if (!bellmanFord(g, s, dist, parentNode, parentArc))
        return false; // negative cycle

    const int INF = numeric_limits<int>::max() / 4;
    int bestDist = INF;
    t = -1; // initialisation : non puit trouve <- -1

    // trouver le puit t avec la plus petite distance parmi les noeuds avec supply negative
    for (int v = 0; v < g.numNodes; ++v)
    {
        // considerer uniquement les noeuds avec supply negative
        if (g.nodes[v].supply < 0 && dist[v] < bestDist)
        {
            // màj du meilleur puit trouve
            bestDist = dist[v];
            t = v;
        }
    }
    return (t != -1);
}

int trouverDeltaMin(const Graph &g,
                    int s, int t,
                    const vector<int> &parentNode,
                    const vector<int> &parentArc)
{
    const int INF = numeric_limits<int>::max() / 4;
    int delta = INF;

    // trouver le minimum capacity residu sur le chemin Wi
    for (int v = t; v != s; v = parentNode[v])
    {
        // remonter le chemin de t à s en utilisant parentNode et parentArc
        // parentNode[v] : noeud prédécesseur de v sur le chemin
        // parentArc[v] : arc utilisé pour atteindre v depuis parentNode[v]
        // donc l'arc (u -> v) est g.adj[u][idx] avec u = parentNode[v] et idx = parentArc[v]
        // on fait un backtracking du puit t vers la source s
        int u = parentNode[v];
        int idx = parentArc[v];
        delta = min(delta, g.adj[u][idx].capacityResidu);
    }
    // considerer les supply des noeuds s et t
    delta = min(delta, g.nodes[s].supply);
    delta = min(delta, -g.nodes[t].supply);
    return delta;
}

void augmenterFlotEtMettreAJourGraphe(Graph &g,
                                      int s, int t,
                                      int delta,
                                      const vector<int> &parentNode,
                                      const vector<int> &parentArc,
                                      int &totalCost)
{
    // augmenter le flot de delta sur le chemin Wi
    for (int v = t; v != s; v = parentNode[v])
    {
        // on fait un backtracking du puit t vers la source s (pareil que dans trouverDeltaMin)
        int u = parentNode[v];
        int idx = parentArc[v];
        Arc &arc = g.adj[u][idx];
        Arc &rev = g.adj[arc.toNode][arc.reverseArcId];

        arc.capacityResidu -= delta;
        rev.capacityResidu += delta;
        totalCost += delta * arc.cost; // màj du coût total
    }
    // màj des supply des noeuds s et t
    g.nodes[s].supply -= delta;
    g.nodes[t].supply += delta;
}

int minCostFlow(Graph &g)
{
    /* On implemente l'algo décrit dans le slide 20 du cours.
       * On cherche à envoyer du flot depuis les noeuds avec supply positive
        vers les noeuds avec supply negative en minimisant le coût total.
       * On utilise Bellman-Ford pour trouver les chemins de coût minimum.
       * On répète jusqu'à ce qu'on ne puisse plus envoyer de flot.
       * On retourne le coût total du flot envoyé.
       * On suppose que le flot de coût minimum est réalisable.
       * On ne gère pas le cas où il n'y a pas de solution réalisable.
    */

    // initialisation
    const int n = g.numNodes;
    int totalCost = 0;
    int iterations = 0;
    const int MAX_ITER = 5000; // pour éviter boucle infinie :) (limite arbitraire)

    while (iterations++ < MAX_ITER)
    {
        bool progress = false; // un checker pour voir si on a fait des augmentations de flot cette itération (utile pour arrêter si plus rien à faire !!!)

        for (int s = 0; s < n; ++s)
        {
            // considérer uniquement les noeuds avec supply positive comme sources car on cherche à envoyer du flot depuis ces noeuds
            if (g.nodes[s].supply <= 0)
                continue;

            vector<int> dist, parentNode, parentArc;                              // pour Bellman-Ford
            int t;                                                                // puit trouvé
            bool okPath = trouverCheminMin(g, s, dist, parentNode, parentArc, t); // trouver le chemin de coût minimum depuis s vers un puit t
            if (!okPath)
                continue; // on n'a pas un cycle négatif ou pas de puit trouvable => passer à la source suivante

            int delta = trouverDeltaMin(g, s, t, parentNode, parentArc);
            if (delta <= 0)
                continue; // pas d'augmentation possible => passer à la source suivante

            augmenterFlotEtMettreAJourGraphe(g, s, t, delta,
                                             parentNode, parentArc,
                                             totalCost);

            progress = true;
        }

        if (!progress) // on a pas pu faire d'augmentation de flot cette itération => on a fini
            break;

        bool done = true; // checker si tous les supply sont nuls (terminaison);
                          // cette étape est optionnelle car on peut aussi se baser sur le fait qu'on ne peut plus faire d'augmentation de flot
        for (const Node &nd : g.nodes)
        {
            if (nd.supply != 0) // il reste du supply à traiter qui sont pas encore à 0
            {
                done = false;
                break;
            }
        }
        if (done)
            break;
    }
    return totalCost;
}

/***********************************************************************/
/******************************** TESTs ********************************/
/***********************************************************************/
bool isFlowValid(const Graph &g)
{
    const int n = g.numNodes;

    // 1. Capacity + reverse consistency
    for (int u = 0; u < n; ++u)
    {
        for (const Arc &arc : g.adj[u])
        {
            int v = arc.toNode;
            int flow = arc.MAXCAP - arc.capacityResidu;

            // Enforce capacity bounds only on forward arcs
            if (arc.MAXCAP > 0)
            {
                if (flow < arc.MINCAP || flow > arc.MAXCAP)
                    return false;
            }

            if (!isValidIndex(g, v) ||
                arc.reverseArcId < 0 ||
                arc.reverseArcId >= static_cast<int>(g.adj[v].size()))
            {
                return false;
            }

            const Arc &rev = g.adj[v][arc.reverseArcId];
            int rev_flow = rev.MAXCAP - rev.capacityResidu;
            if (flow + rev_flow != 0)
                return false;
        }
    }

    // 2. Flow conservation
    vector<int> balance(n, 0);
    for (int u = 0; u < n; ++u)
    {
        for (const Arc &arc : g.adj[u])
        {
            if (arc.MAXCAP == 0)
                continue; // Skip reverse arcs
            int v = arc.toNode;
            int flow = arc.MAXCAP - arc.capacityResidu;
            balance[u] -= flow;
            balance[v] += flow;
        }
    }

    // 3. Check against original supplies (outflow - inflow = -supply)
    for (int i = 0; i < n; ++i)
    {
        if (balance[i] != -g.nodes[i].originalSupply)
        {
            cout << "JUST A PRINT: Balance fail at node " << i
                 << ": balance=" << balance[i]
                 << " vs -originalSupply=" << -g.nodes[i].originalSupply << endl;
            return false;
        }
    }
    return true;
}

// Define the macro for assert (code used from https://www.geeksforgeeks.org/cpp/how-to-add-message-to-assert-in-cpp/)
#define ASSERT(condition, message)     \
    do                                 \
    {                                  \
        assert(condition && #message); \
    } while (0)

void checkCoutMin()
{
    string filename = "./instances/example.dat";
    Graph g = readGraph(filename);
    ASSERT(minCostFlow(g) == 209195 && isFlowValid(g) == true, "Failed on example.dat");

    filename = "./instances/example_small_1.dat";
    g = readGraph(filename);
    ASSERT(minCostFlow(g) == 141486 && isFlowValid(g) == true, "Failed on example_small_1.dat");
    filename = "./instances/example_small_2.dat";
    g = readGraph(filename);
    ASSERT(minCostFlow(g) == 595779 && isFlowValid(g) == true, "Failed on example_small_2.dat");
    filename = "./instances/example_small_3.dat";
    g = readGraph(filename);
    ASSERT(isFlowValid(g) == false, "Failed on example_small_3.dat"); // NOT POSSIBLE
    filename = "./instances/example_small_4.dat";
    g = readGraph(filename);
    ASSERT(minCostFlow(g) == 9108 && isFlowValid(g) == true, "Failed on example_small_4.dat");

    filename = "./instances/example_medium_1.dat";
    g = readGraph(filename);
    ASSERT(minCostFlow(g) == 332662 && isFlowValid(g) == true, "Failed on example_medium_1.dat");
    filename = "./instances/example_medium_2.dat";
    g = readGraph(filename);
    ASSERT(minCostFlow(g) == 3320788 && isFlowValid(g) == true, "Failed on example_medium_2.dat");
    filename = "./instances/example_medium_3.dat";
    g = readGraph(filename);
    ASSERT(minCostFlow(g) == 313037611 && isFlowValid(g) == true, "Failed on example_medium_3.dat");
    filename = "./instances/example_medium_4.dat";
    g = readGraph(filename);
    ASSERT(minCostFlow(g) == 6517946 && isFlowValid(g) == true, "Failed on example_medium_4.dat");
    filename = "./instances/example_medium_5.dat";
    g = readGraph(filename);
    ASSERT(minCostFlow(g) == 2157318 && isFlowValid(g) == true, "Failed on example_medium_5.dat");

    filename = "./instances/example_large_1.dat";
    g = readGraph(filename);
    ASSERT(minCostFlow(g) == 919960 && isFlowValid(g) == true, "Failed on example_large_1.dat");
    filename = "./instances/example_large_2.dat";
    g = readGraph(filename);
    ASSERT(minCostFlow(g) == 9085629 && isFlowValid(g) == true, "Failed on example_large_2.dat");

    // filename = "./instances/example_verylarge_1.dat";
    // g = readGraph(filename);
    // ASSERT(minCostFlow(g) == && isFlowValid(g) == true, "Failed on example_verylarge_1.dat");
    // filename = "./instances/example_verylarge_2.dat";
    // g = readGraph(filename);
    // ASSERT(minCostFlow(g) == && isFlowValid(g) == true, "Failed on example_verylarge_2.dat");

    cout << "All cout min tests passed! {^_^}" << endl;
}

/***********************************************************************/
/************************ PRINT SOLUTION *******************************/
/***********************************************************************/

void printSolution(string filename)
{
    /* Prints "Coût optimal"*/
    Graph g = readGraph(filename);

    cout << "Initial graph loaded" << endl;

    int totalCost = minCostFlow(g);

    cout << "Flow valid? " << (isFlowValid(g) ? "YES" : "NO") << endl;
    cout << "Coût optimal: " << totalCost << endl;

    // cout << "Flots sur chaque arc:" << endl;
    // for (int u = 0; u < g.numNodes; ++u)
    // {
    //     for (const Arc &arc : g.adj[u])
    //     {
    //         if (arc.MAXCAP > 0)
    //         { // Only forward arcs (not reverse)
    //             int flow = arc.MAXCAP - arc.capacityResidu;
    //             cout << "Arc " << (u + 1) << " -> " << (arc.toNode + 1)
    //                  << ": flot = " << flow << endl;
    //         }
    //     }
    // }
}

/***********************************************************************/
/******************************** MAIN *********************************/
/***********************************************************************/

int main()
{
    cout << "HELLO :), FLOT_COUT_MIN" << endl;

    checkCoutMin();

    try // Important for catching file read errors implemented in the above functions
    {
        // string filename = "./example0Min.txt";
        // string filename = "./instances/example.dat";
        // string filename = "./instances/example_small_3.dat";

        // printSolution(filename);

        // Graph g = readGraph(filename);

        // vector<int> dist, parentNode, parentArc;

        // // Initialize distances to 0 for all nodes to check for negative cycles everywhere
        // const int n = g.numNodes;
        // const int INF = numeric_limits<int>::max() / 4;
        // dist.assign(n, 0); // Start all distances at 0
        // parentNode.assign(n, -1);
        // parentArc.assign(n, -1);

        // // Relax edges n times (consider ALL arcs including reverse arcs)
        // for (int k = 0; k < n; ++k)
        // {
        //     for (int u = 0; u < n; ++u)
        //     {
        //         for (int idx = 0; idx < static_cast<int>(g.adj[u].size()); ++idx)
        //         {
        //             const Arc &arc = g.adj[u][idx];
        //             // Don't skip any edges - check all including reverse arcs

        //             int v = arc.toNode;
        //             if (dist[u] + arc.cost < dist[v])
        //             {
        //                 dist[v] = dist[u] + arc.cost;
        //                 parentNode[v] = u;
        //                 parentArc[v] = idx;
        //             }
        //         }
        //     }
        // }

        // // Check if we can still relax any edge (negative cycle detection)
        // bool hasNegCycle = false;
        // for (int u = 0; u < n && !hasNegCycle; ++u)
        // {
        //     for (int idx = 0; idx < static_cast<int>(g.adj[u].size()); ++idx)
        //     {
        //         const Arc &arc = g.adj[u][idx];
        //         // Check all edges including reverse arcs

        //         int v = arc.toNode;
        //         if (dist[u] + arc.cost < dist[v])
        //         {
        //             hasNegCycle = true;
        //             break;
        //         }
        //     }
        // }
        // cout << "Existence de cycle négatif ? " << (hasNegCycle ? "OUI" : "NON") << endl;
    }
    catch (const exception &e)
    {
        cerr << "ERROR: " << e.what() << endl;
        return 1;
    }

    return 0;
}
