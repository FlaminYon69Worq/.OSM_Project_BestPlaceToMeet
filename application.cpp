// application.cpp <Starter Code>
// <David Mendoza>
//
// University of Illinois at Chicago
// CS 251: Fall 2021
// Project #7 - Openstreet Maps
//
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//
#include <algorithm>
#include <iostream>
#include <iomanip>  /*setprecision*/
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <limits>
#include <stack>
#include "tinyxml2.h"
#include "dist.h"
#include "osm.h"
#include "graph.h"

using namespace std;
using namespace tinyxml2;

//
// Implement your creative component application here
// TO DO: add arguments
//
const double INF = numeric_limits<double>::max();
class prioritize{
  public:
      bool operator() (pair<long long, double> p,
      pair<long long, double > p1) const {
        return p.second > p1.second;
      }
};

void creative() {
    cout << "undone :(" << endl;
}

// searchBuilding: function to find building in vector
BuildingInfo searchBuilding(string query, vector<BuildingInfo>& Buildings) {
  BuildingInfo found;
  for (auto b : Buildings) {
  // loop thru buildings
    if (b.Abbrev.find(query)
             != string::npos || b.Fullname.find(query)
             != string::npos) {
      found = b;  // word found, set found as the cur building
      return found;
    }
  }
  return found;
}

// nearestBuilding: takes center, finds min distance between two pts
BuildingInfo nearestBuilding(BuildingInfo& t1, BuildingInfo& t2, vector<BuildingInfo>& Buildings) {
  BuildingInfo nearest;
  double min = INF;  // minimum distance
  double dist;  // distances we are calculating
  Coordinates midpoint =  // down below
  centerBetween2Points(t1.Coords.Lat, t1.Coords.Lon,
                               t2.Coords.Lat, t2.Coords.Lon);
  // loop thru buildings, min algorithm
  for (auto b : Buildings) {
    dist = distBetween2Points(midpoint.Lat, midpoint.Lon,
                             b.Coords.Lat,
                             b.Coords.Lon);
    if (dist < min) {
      min = dist;
      nearest = b;
    }
  }
  return nearest;
}

// nearestBuildingwSet: takes center, finds min distance between two pts
BuildingInfo nearestBuildingwSet(BuildingInfo& t1, BuildingInfo& t2, vector<BuildingInfo>& Buildings, set<string> s) {
  BuildingInfo nearest;
  double min = INF;  // minimum distance
  double dist;  // distances we are calculating
  Coordinates midpoint =  // down below
  centerBetween2Points(t1.Coords.Lat, t1.Coords.Lon,
                               t2.Coords.Lat, t2.Coords.Lon);
  // loop thru buildings
  for (auto b : Buildings) {
    dist = distBetween2Points(midpoint.Lat, midpoint.Lon,
                             b.Coords.Lat,
                             b.Coords.Lon);
    if (dist < min && s.count(b.Fullname) == 0) {
      min = dist;
      nearest = b;
    }
  }
  return nearest;
}

// // nearestBuilding: uses building to find nearest node, returns ID
long long nearestNode(BuildingInfo& b, map<long long, Coordinates>& Nodes,
vector<FootwayInfo>& Footways) {
  double dist, min = INF;
  long long returnID;
  for (auto& f : Footways) {
    for (unsigned int i = 0; i < f.Nodes.size(); i++) {  // Nodes vector
      dist = distBetween2Points(Nodes[f.Nodes[i]].Lat,
                                      Nodes[f.Nodes[i]].Lon,
                                      b.Coords.Lat, b.Coords.Lon);
      if (dist < min) {  // minimum algorithm
        min = dist;
        returnID = Nodes[f.Nodes[i]].ID;
      }
    }
  }
  return returnID;
}

// popUnvisited: while loop of djikstra algorithm
// helper fxn
void popUnvisited(
  priority_queue<pair<long long, double>,
  vector<pair<long long, double>>, prioritize> pq,
  graph<long long, double>& G, vector<long long> visited, set<long long> discovered,
  long long endV, map<long long, long long>& pred,
  map<long long, double>& distances) {
    while (!pq.empty()) {
     // Visit vertex with minimum distance from startV
     pair<long long, double> w;
     w.first = pq.top().first;
     w.second = pq.top().second;
     pq.pop();
    // logic:
    //  if (currentV->distance == infinity)
    //   break;
    // else if (currentV has been visited)
    //   continue;
    // else
    //   visit currentV
    //
     if (w.second == INF) {
         break;
     } else if (w.first == endV) {
         visited.push_back(endV);
         break;
     } else if (discovered.count(w.first) == 1) {
         continue;
     } else {  // add to storage of discovered vertices
         discovered.insert(w.first);
         visited.push_back(w.first);
     }
     set<long long> Neighbors = G.neighbors(w.first);
     for (auto&j : Neighbors) {
       double weight;
       G.getWeight(w.first, j, weight);
       double alt = distances[w.first] + weight;
       // If shorter path from startV to adjV is found,
       // update adjV's distance and predecessor
       if (alt < distances[j]) {
         distances[j] = alt;
         pred[j] = w.first;
         pair<long long, double> l;
         l.first = j;
         l.second = alt;
         pq.push(l);
      }
    }
  }
}

// djikstra's algorithm
void DjikstraAlg(graph<long long, double>& G, long long startV, long long endV,
map<long long, double>& distances, map<long long, long long>& pred) {
  vector<long long> visited;
  set<long long> discovered;
  priority_queue<pair<long long, double>,
  vector<pair<long long, double>>, prioritize> pq;
  vector<long long> allVerts = G.getVertices();
  //
  //  for each vertex currentV in graph {
  //     currentV⇢distance = Infinity
  //     currentV⇢predV = 0
  //     Enqueue currentV in unvisitedQueue
  //  }
  for (auto& v : allVerts) {
    pair<long long, double> p;
    p.first = v;
    p.second = INF;
    pq.push(p);
    distances[v] = INF;
    pred[v] = 0;
  }
  // startV has a distance of 0 from itself
  pair<long long, double> firstPair;
  firstPair.first = startV;
  firstPair.second = 0;
  pq.push(firstPair);
  distances[startV] = 0;
  popUnvisited(pq, G, visited, discovered, endV, pred, distances);
}

// getting the path of a person
vector<long long> ms11Path(graph<long long, double>& G, long long startV,
long long endV, map<long long, double>& distances, map<long long, long long>& pred) {
  stack<long long> s;
  vector<long long> path;
  long long vert = endV;
  while (vert != 0) {
    s.push(vert);
    vert = pred[vert];
  }
  while (!s.empty()) {
    long long node = s.top();  // current node
    path.push_back(node);  // push into the path
    s.pop();  // node visited, pop it
  }
  return path;
}

// prints person 1, 2 and the destinations info
void printLocations(BuildingInfo& b1, BuildingInfo& b2, BuildingInfo& NB,
map<long long, Coordinates>& Nodes, long long n1, long long n2,
long long dest) {
  cout << "Person 1's point:" << endl;
  cout << " " << b1.Fullname << endl;
  cout <<" (" << b1.Coords.Lat << ", " << b1.Coords.Lon<< ")" <<endl;
  cout << "Person 2's point:"<< endl;
  cout << " " << b2.Fullname << endl;
  cout << " (" << b2.Coords.Lat << ", " << b2.Coords.Lon << ")" << endl;
  cout << "Destination Building:" << endl;
  cout << " " << NB.Fullname << endl;
  cout << " (" << NB.Coords.Lat << ", " << NB.Coords.Lon << ")" << endl;
  cout << endl;
  cout << "Nearest P1 node:" << endl;
  cout << " " << n1 << endl;
  cout << " (" << Nodes[n1].Lat << ", " << Nodes[n1].Lon << ")" << endl;
  cout << "Nearest P2 node:" << endl;
  cout << " " << n2 << endl;
  cout << " (" << Nodes[n2].Lat << ", " << Nodes[n2].Lon<< ")" <<endl;
  cout << "Nearest destination node:" << endl;
  cout << dest << endl;
  cout << " (" << Nodes[dest].Lat << ", "<< Nodes[dest].Lon << ")" << endl;
  cout << endl;
}

// application
void application(graph<long long, double>& G,
    map<long long, Coordinates>& Nodes, vector<FootwayInfo>& Footways,
    vector<BuildingInfo>& Buildings) {
  string person1Building, person2Building;

  cout << endl;
  cout << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, person1Building);

  while (person1Building != "#") {
    cout << "Enter person 2's building (partial name or abbreviation)> ";
    getline(cin, person2Building);
    // search for person1 and 2 building, and midpoint, if possible
    BuildingInfo b1 = searchBuilding(person1Building, Buildings);
    BuildingInfo b2 = searchBuilding(person2Building, Buildings);
    BuildingInfo NB = nearestBuilding(b1, b2, Buildings);
    // search nearest node of buildings
    long long n1 = nearestNode(b1, Nodes, Footways);
    long long n2 = nearestNode(b2, Nodes, Footways);
    long long dest = nearestNode(NB, Nodes, Footways);
    // declare paths of person 1, 2 and btwn 1 and 2 respectively
    // 1 to dest, 2 to dest, and 3 is btwn 1 and 2
    vector<long long> p1;
    vector<long long> p2;
    vector<long long> p3;
    // predecessor and distances of person 1 and 2
    map<long long, long long> pred1;
    map<long long, long long> pred2;
    map<long long, double> d1;
    map<long long, double> d2;
    // compute djikstra and get path
    DjikstraAlg(G, n1, dest , d1, pred1);
    p1 = ms11Path(G, n1, dest , d1, pred1);
    // do the same with the person 2
    DjikstraAlg(G, n2, dest , d2, pred2);
    p2 = ms11Path(G, n2, dest, d2, pred2);
    // checks if n1 to n2 is reachable, 3rd djikstra and path
    map<long long, long long> pred3;
    map<long long, double> d3;
    DjikstraAlg(G, n1, n2, d3, pred3);
    p3 = ms11Path(G, n1, n2, d3, pred3);

    set<string> buildings;  // possible destinations
    // 0.0 is not a valid coordinate, unless there is a place
    // in the North Pole lol
    if (b1.Coords.Lat == 0.0 && b1.Coords.Lon == 0.0) {
        cout << "Person 1's building not found" << endl;
    } else if (b2.Coords.Lat == 0.0 && b2.Coords.Lon == 0.0) {
        cout << "Person 2's building not found" << endl;
    } else {  // valid coordinates
      // print out everyone's location, and node
      // also print out destination location and node
      printLocations(b1, b2, NB, Nodes, n1, n2, dest);
      if (d3[n2] >= INF) {  // path from n1 to n2 unreachable
        cout << "Sorry, destination unreachable." << endl;
        cout << endl;
      } else if (d1[dest] >= INF || d2[dest] >= INF) {  // down below
          cout << "At least one person was unable to reach the destination building. Finding next closest building..." << endl;
          cout << endl;
          buildings.insert(NB.Fullname);  // possible destination marked off
          while (true) {  // while someone cannot reach
              // declare new destination building and node
              NB = nearestBuildingwSet(b1, b2, Buildings, buildings);
              dest = nearestNode(NB, Nodes, Footways);
              if (!(d1[dest] >= INF) || !(d2[dest] >= INF)) {
                break;
              } else {
                  buildings.insert(NB.Fullname);  // destination marked off
                  cout << "New destination building:" << endl;
                  cout << " " << NB.Fullname << endl;
                  cout << " (" << NB.Coords.Lat << ", " << NB.Coords.Lon << ")" << endl;
                  cout << "Nearest destination node:" << endl;
                  cout << dest << endl;
                  cout << " (" << Nodes[dest].Lat << ", "<< Nodes[dest].Lon << ")" << endl;
                  cout << endl;
                  cout << "At least one person was unable to reach the destination building. Finding next closest building..." << endl;
                  cout << endl;
                }
          }  // someone or both can finally reach
          cout << "New destination building:" << endl;
          cout << " " << NB.Fullname << endl;
          cout << " (" << NB.Coords.Lat << ", " << NB.Coords.Lon << ")" << endl;
          cout << "Nearest destination node:" << endl;
          cout << dest << endl;
          cout << " (" << Nodes[dest].Lat << ", "<< Nodes[dest].Lon << ")" << endl;
          // declare new paths and distances to new destination
          vector<long long> new_v;
          vector<long long> new_v2;
          map<long long, double> new_d1;
          map<long long, double> new_d2;
          map<long long, long long> new_pred;
          DjikstraAlg(G, n1, dest, new_d1, new_pred);
          new_v = ms11Path(G, n1, dest, new_d1, new_pred);
          map<long long, long long> new_pred2;
          DjikstraAlg(G, n2, dest, new_d2, new_pred2);
          new_v2 = ms11Path(G, n2, dest, new_d2, new_pred2);
          // print out each person's paths
          cout << "Person 1's distance to dest: " << new_d1.at(dest) << " miles" << endl;
          int i = 0;
          int size = new_v.size();
          cout << "Path: ";
          for (auto x : new_v) {
              if (size-1 == i) {
                break;
              }
              cout << x << "->";
              i++;
          }
          cout << dest << endl;
          cout << endl;
          cout << "Person 2's distance to dest: " << new_d2.at(dest) << " miles" << endl;
          int j = 0;
            size = new_v2.size();
            cout << "Path: ";
          for (auto y : new_v2) {
              if (size-1 == j) {
                break;
              }
              cout << y << "->";
              j++;
          }
          cout <<dest << endl;
      } else {  // destination reachable by both, print out paths
        cout << "Person 1's distance to dest: " << d1.at(dest) << " miles" << endl;
        int i = 0;
        int size = p1.size();
        cout << "Path: ";
        for (auto x : p1) {
            if (size-1 == i) {
              break;
            }
            cout << x << "->";
            i++;
        }
        cout << dest << endl;
        cout << endl;
      
        cout << "Person 2's distance to dest: " << d2.at(dest) << " miles" << endl;
        int j = 0;
          size = p2.size();
          cout << "Path: ";
        for (auto y : p2) {
            if (size-1 == j) {
              break;
            }
            cout << y << "->";
            j++;
        }
        cout << dest << endl;
      }
    }
    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);
  }
}

// helper function in main to make graph, save space
void makeGraph(graph<long long, double>& G, map<long long, Coordinates>& Nodes,
vector<FootwayInfo>& Footways) {
  for (auto e : Nodes) {
    G.addVertex(e.first);
}

  for (auto& l : Footways) {  // add edges
      for (unsigned int j = 0; j < l.Nodes.size()-1; j++) {
        double distance = distBetween2Points(Nodes[l.Nodes.at(j)].Lat,
                        Nodes[l.Nodes.at(j)].Lon, Nodes[l.Nodes.at(j+1)].Lat,
                        Nodes[l.Nodes.at(j+1)].Lon);
        G.addEdge(l.Nodes.at(j), l.Nodes.at(j+1), distance);
        G.addEdge(l.Nodes.at(j+1), l.Nodes.at(j), distance);
      }
  }
}

int main() {
  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates>  Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo>          Footways;
  // info about each building, in no particular order
  vector<BuildingInfo>         Buildings;
  XMLDocument                  xmldoc;
  // make graph
  graph<long long, double> G;

  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);


  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }

  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);

  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  //
  // Stats
  //
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;

  makeGraph(G, Nodes, Footways);

  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;

  //
  // Menu
  //
  string userInput;
  cout << "Enter \"a\" for the standard application or "
        << "\"c\" for the creative component application> ";
  getline(cin, userInput);
  if (userInput == "a") {
    application(G, Nodes, Footways, Buildings);
  } else if (userInput == "c") {
    creative();
  }
  //
  // done:
  //
  cout << "** Done **" << endl;
  return 0;
}
