#include<bits/stdc++.h>
#define FOR(i,start,end) for(size_t i = (start); i < (end); i++)
#define pi (2*acos(0.0))
#define MAX 1000
#define eps 0.00001
using namespace std;

struct Point
{
	double x, y;
	int id;

	Point() { id = -2; }
	Point(double a, double b, int i = -1)
	{
		x = a;
		y = b;
		id = i;
	}
};

struct vect
{
	double x, y;

	vect() {}
	vect(Point a, Point b)
	{
		x = b.x - a.x;
		y = b.y - a.y;
	}

	void normalize()
	{
		double mod;
		mod = sqrt(x*x + y*y);
		x = x/mod;
		y = y/mod;
	}

	void rotation(double degree)
	{
		degree = degree * (pi/180); // convert to radian
		double px = x, py = y;
		x = px*cos(degree) - py*sin(degree);
		y = px*sin(degree) + py*cos(degree);
	}

	// (this x v) : ">0" then this is right of v, "<0" this is left of v.
	double crossProduct(vect v)
	{
		return x * (v.y) - y * (v.x);
	}

	double dotProduct(vect v)
	{
		return x * (v.x) + y * (v.y);
	}
};

// ac x ab : ">0" then ac is right of ab, "<0" ac is left of ab.
double crossProduct(Point a, Point b, Point c)
{
	double ac_x = c.x-a.x;
	double ac_y = c.y-a.y;
	double ab_x = b.x-a.x;
	double ab_y = b.y-a.y;
	return ac_x*ab_y - ac_y*ab_x;
}


double dotProduct(vect a, vect b)
{
	return a.x * b.x + a.y * b.y;
}

double angleBetweenVect(vect a, vect b)
{
	a.normalize();
	b.normalize();
	return acos(dotProduct(a, b));
}

inline bool operator == (const Point & p1, const Point & p2)
{
	return abs(p1.x - p2.x) < eps && abs(p1.y - p2.y) < eps;
}

struct DelaunayEdge
{
	Point startP, endP, uP, dP;
	int id;
	bool isValid;

	DelaunayEdge() {}
	DelaunayEdge(Point start = Point(), Point end = Point(), Point up = Point(), Point down = Point())
	{
		startP = start;
		endP = end;
		uP = up;
		dP = down;
	}
};

ostream & operator << (ostream & stream, const Point& p)
{
	stream << p.x << " " << p.y;
	return stream;
}

ostream & operator << (ostream & stream, const DelaunayEdge & D)
{
	stream << "start : " << D.startP << " end : " << D.endP << " up : " << D.uP << " down : " << D.dP << " edge id " << D.id<< endl;
	return stream;
}

// variables
int numberOfPoints;
int numberOfInputEdges = 3;
vector<Point> givenPoints;
vector<DelaunayEdge> edgeList;
vector< vector <int> > mapping;
ifstream inputFile;
ofstream outputFile;

// functions
void initialize();
bool isInTriangle(Point & p, Point & a, Point & b, Point & c);
bool isOnEdge(Point & p, Point & a, Point & b);
void legalizeEdge(int pointId, int edgeId);

int main()
{
	ios_base::sync_with_stdio(0);
	cin.tie(0);
	cout.tie(0);

	initialize();
	int invalidId_1 = givenPoints[givenPoints.size()-2].id;
	int invalidId_2 = givenPoints[givenPoints.size()-1].id;

	// main procedure
	FOR(i, 1, givenPoints.size()-2)
	{
		Point Pr = givenPoints[i];
		FOR(j, 0, edgeList.size())
		{
			if (edgeList[j].isValid)
			{
				Point S = edgeList[j].startP;
				Point E = edgeList[j].endP;
				Point U = edgeList[j].uP;
				Point D = edgeList[j].dP;

				if (isInTriangle(Pr, S, E, U))
				{
					// in the upper triangle
					// creating new 3 edges
					// upper left
					DelaunayEdge edge1(Pr, S, U, E);
					edge1.isValid = 1;
					edge1.id = numberOfInputEdges++;
					edgeList.push_back(edge1);
					mapping[edge1.startP.id][edge1.endP.id] = edge1.id;
					mapping[edge1.endP.id][edge1.startP.id] = edge1.id;

					// upper right
					DelaunayEdge edge2(Pr, E, U, S);
					edge2.isValid = 1;
					edge2.id = numberOfInputEdges++;
					edgeList.push_back(edge2);
					mapping[edge2.startP.id][edge2.endP.id] = edge2.id;
					mapping[edge2.endP.id][edge2.startP.id] = edge2.id;

					// up
					DelaunayEdge edge3(Pr, U, S, E);
					edge3.isValid = 1;
					edge3.id = numberOfInputEdges++;
					edgeList.push_back(edge3);
					mapping[edge3.startP.id][edge3.endP.id] = edge3.id;
					mapping[edge3.endP.id][edge3.startP.id] = edge3.id;

					// find edge 1,2,3
					// edge 1
					int edgeSEid = mapping[S.id][E.id];//edgeList[j].id;
					// edge 2
					int edgeEUid = mapping[E.id][U.id];
					// edge 3
					int edgeSUid = mapping[S.id][U.id];

					if (edgeSEid == -1 || edgeEUid == -1 || edgeSUid == -1)
						cout << "Error choosing edge at legalize\n";

					// update old edges
					// for edge 1
					edgeList[edgeSEid].uP = Pr;
					// for edge 2
					if (edgeList[edgeSEid].startP.id == edgeList[edgeEUid].uP.id)
						edgeList[edgeEUid].uP = Pr;
					else if (edgeList[edgeSEid].startP.id == edgeList[edgeEUid].dP.id)
						edgeList[edgeEUid].dP = Pr;
					// for edge 3
					if (edgeList[edgeSEid].endP.id == edgeList[edgeSUid].uP.id)
						edgeList[edgeSUid].uP = Pr;
					else if (edgeList[edgeSEid].endP.id == edgeList[edgeSUid].dP.id)
						edgeList[edgeSUid].dP = Pr;

					// calling legalize
					legalizeEdge(Pr.id, edgeSEid);
					legalizeEdge(Pr.id, edgeEUid);
					legalizeEdge(Pr.id, edgeSUid);

					break;
				}
				else if (isInTriangle(Pr, S, E, D))
				{
					// in the lower triangle
					// creating new 3 edges
					// lower left
					DelaunayEdge edge1(Pr, S, E, D);
					edge1.isValid = 1;
					edge1.id = numberOfInputEdges++;
					edgeList.push_back(edge1);
					mapping[edge1.startP.id][edge1.endP.id] = edge1.id;
					mapping[edge1.endP.id][edge1.startP.id] = edge1.id;

					// lower right
					DelaunayEdge edge2(Pr, E, D, S);
					edge2.isValid = 1;
					edge2.id = numberOfInputEdges++;
					edgeList.push_back(edge2);
					mapping[edge2.startP.id][edge2.endP.id] = edge2.id;
					mapping[edge2.endP.id][edge2.startP.id] = edge2.id;

					// low
					DelaunayEdge edge3(Pr, D, S, E);
					edge3.isValid = 1;
					edge3.id = numberOfInputEdges++;
					edgeList.push_back(edge3);
					mapping[edge3.startP.id][edge3.endP.id] = edge3.id;
					mapping[edge3.endP.id][edge3.startP.id] = edge3.id;

					// find edge 1,2,3
					// edge 1
					int edgeSEid = edgeList[j].id;
					// edge 2
					int edgeEDid = mapping[E.id][D.id];
					// edge 3
					int edgeSDid = mapping[S.id][D.id];

					if (edgeSEid == -1 || edgeEDid == -1 || edgeSDid == -1)
						cout << "Error choosing edge at legalize\n";

					// update old edges
					// for edge 1
					edgeList[edgeSEid].dP = Pr;
					// for edge 2
					if (edgeList[edgeSEid].startP.id == edgeList[edgeEDid].uP.id)
						edgeList[edgeEDid].uP = Pr;
					else if (edgeList[edgeSEid].startP.id == edgeList[edgeEDid].dP.id)
						edgeList[edgeEDid].dP =  Pr;
					// for edge 3
					if (edgeList[edgeSEid].endP.id == edgeList[edgeSDid].uP.id)
						edgeList[edgeSDid].uP = Pr;
					else if (edgeList[edgeSEid].endP.id == edgeList[edgeSDid].dP.id)
						edgeList[edgeSDid].dP = Pr;

					// calling legalize
					legalizeEdge(Pr.id, edgeSEid);
					legalizeEdge(Pr.id, edgeEDid);
					legalizeEdge(Pr.id, edgeSDid);

					break;
				}
				else if (isOnEdge(Pr, S, E))
				{
					// on the edge
					// remove current edge
					edgeList[j].isValid = 0;
                    mapping[edgeList[j].startP.id][edgeList[j].endP.id] = -1;
                    mapping[edgeList[j].endP.id][edgeList[j].startP.id] = -1;

                    // create new edge
                    DelaunayEdge edge1(Pr, S, U, D);
                    edge1.isValid = 1;
                    edge1.id = numberOfInputEdges++;
                    edgeList.push_back(edge1);
                    mapping[edge1.startP.id][edge1.endP.id] = edge1.id;
					mapping[edge1.endP.id][edge1.startP.id] = edge1.id;

                    DelaunayEdge edge2(Pr, D, E, S);
					edge2.isValid = 1;
					edge2.id = numberOfInputEdges++;
					edgeList.push_back(edge2);
					mapping[edge2.startP.id][edge2.endP.id] = edge2.id;
					mapping[edge2.endP.id][edge2.startP.id] = edge2.id;

					DelaunayEdge edge3(Pr, E, U, D);
					edge3.isValid = 1;
					edge3.id = numberOfInputEdges++;
					edgeList.push_back(edge3);
					mapping[edge3.startP.id][edge3.endP.id] = edge3.id;
					mapping[edge3.endP.id][edge3.startP.id] = edge3.id;

					DelaunayEdge edge4(Pr, U, S, E);
					edge4.isValid = 1;
					edge4.id = numberOfInputEdges++;
					edgeList.push_back(edge4);
					mapping[edge4.startP.id][edge4.endP.id] = edge4.id;
					mapping[edge4.endP.id][edge4.startP.id] = edge4.id;

					// find edge 1,2,3,4
                    int edgeSDid = mapping[S.id][D.id];
                    int edgeEDid = mapping[D.id][E.id];
                    int edgeEUid = mapping[E.id][U.id];
                    int edgeSUid = mapping[S.id][U.id];

                    if (edgeSDid == -1 || edgeEDid == -1 || edgeEUid == -1 || edgeSUid == -1)
						cout << "Error choosing edge at legalize\n";

					legalizeEdge(Pr.id, edgeSDid);
					legalizeEdge(Pr.id, edgeEDid);
					legalizeEdge(Pr.id, edgeEUid);
					legalizeEdge(Pr.id, edgeSUid);

					break;
				}
				else{
					//cout << "error at " << givenPoints[i] << " not on edge or triangle" << endl;
					//break;
				}
			}

		}
	}

	// removing the added invalid edges
	FOR(i, 0, edgeList.size())
	{
		if ((edgeList[i].startP.id == -1 || edgeList[i].startP.id == -2) && (edgeList[i].endP.id == -1 || edgeList[i].endP.id == -2)) {
			edgeList[i].isValid = false;
		}
		if ((edgeList[i].startP.id == invalidId_1 || edgeList[i].startP.id == invalidId_2) || (edgeList[i].endP.id == invalidId_1 || edgeList[i].endP.id == invalidId_2)) {
			edgeList[i].isValid = false;
		}
	}

	FOR(i, 0, edgeList.size())
	{
		if (edgeList[i].isValid) {
			outputFile << edgeList[i].startP << " " << edgeList[i].endP;
            if (i != edgeList.size()-1) outputFile << "\n";
		}
	}

	return 0;
}

void initialize()
{
	freopen("output.txt", "w", stdout);
	inputFile.open("input.txt");
	outputFile.open("inputDelaunay.txt");
	inputFile >> numberOfPoints;
	givenPoints.assign(numberOfPoints, Point());
	FOR(i, 0, numberOfPoints)
	{
		inputFile >> givenPoints[i].x >> givenPoints[i].y;
	}

	// sort according to y and x basis
	sort(givenPoints.begin(), givenPoints.end(),[](const Point& upPoint, const Point& downPoint)
	{
		return upPoint.y == downPoint.y ? upPoint.x > downPoint.x : upPoint.y > downPoint.y ;
	});

	// p0 is the topmost rightmost point
	Point p0 = *(givenPoints.begin());
	random_shuffle(givenPoints.begin() + 1, givenPoints.end());

	// find rightmost
	Point pR = givenPoints[1];
	FOR(i, 2, givenPoints.size())
	{
		//  the chosen pR is not at right
		if (crossProduct(p0, givenPoints[i], pR) < 0) {
			pR = givenPoints[i];
		}
	}
	// find leftmost
	Point pL = givenPoints[1];
	FOR(i, 2, givenPoints.size())
	{
		//  the chosen pL is not at left
		if (crossProduct(p0, pL, givenPoints[i]) < 0) {
			pL = givenPoints[i];
		}
	}

	vect leftSide(p0, pL);
	vect rightSide(p0, pR);
	leftSide.normalize();
	rightSide.normalize();
	leftSide.rotation(10);
	rightSide.rotation(-10);

	// leftmost pl_2, righmost pr_1
	Point pr_1(int(p0.x + rightSide.x * MAX), int(p0.y + rightSide.y * MAX));
	Point pl_2(int(p0.x + leftSide.x * MAX), int(p0.y + leftSide.y * MAX));

	givenPoints.push_back(pr_1);
	givenPoints.push_back(pl_2);
	FOR(i, 0, givenPoints.size())
	{
		givenPoints[i].id = i;
	}

	FOR(i, 0, givenPoints.size())
	{
		cout << givenPoints[i] << endl;
	}
	cout << endl;

	// need to check mappings, error occuring or not
	mapping.assign(givenPoints.size(), vector<int>());
	FOR(i, 0, givenPoints.size())
	{
		mapping[i].assign(givenPoints.size(), -1);
	}

	// pushing values in edgeList
	p0 = givenPoints[0];
	pr_1 = givenPoints[givenPoints.size()-2];
	pl_2 = givenPoints[givenPoints.size()-1];
	edgeList.push_back(DelaunayEdge(pr_1, pl_2, p0));
	edgeList.push_back(DelaunayEdge(p0, pr_1, pl_2));
	edgeList.push_back(DelaunayEdge(p0, pl_2, pr_1));
	// mapping in the ara
	FOR(i, 0, edgeList.size())
	{
		edgeList[i].id = i;
		edgeList[i].isValid = 1;
		mapping[edgeList[i].startP.id][edgeList[i].endP.id] = edgeList[i].id; // edgeList[i].id;
		mapping[edgeList[i].endP.id][edgeList[i].startP.id] = edgeList[i].id;
	}
}

bool isInTriangle(Point & p, Point & a, Point & b, Point & c)
{

	if (a.id == -2 || b.id == -2 || c.id == -2 || p.id == -2) {
		//cout << "error\na " << a << " b " << b << " c " << c << " p " << p << endl;
		return 0;
	}
	vect ab(a, b), bc(b, c), ca(c, a);
	vect ap(a, p), bp(b, p), cp(c, p);

	double abc, abp, bcp, cap;

	abp = 0.5 * abs(ab.crossProduct(ap));
	bcp = 0.5 * abs(bc.crossProduct(bp));
	cap = 0.5 * abs(ca.crossProduct(cp));
	abc = 0.5 * abs(-ab.crossProduct(bc));

	return abs(abc - (abp + bcp + cap)) < eps;
}

bool isOnEdge(Point & p, Point & a, Point & b)
{
	return crossProduct(a, p, b) == 0;  // ab x ap == 0
}

void legalizeEdge(int pointId, int edgeId)
{
	if (!edgeList[edgeId].isValid) return;

	Point S = edgeList[edgeId].startP;
	Point E = edgeList[edgeId].endP;
	Point U = edgeList[edgeId].uP;
	Point D = edgeList[edgeId].dP;

	// if no D or U then valid
	if (D.id == -2 || U.id == -2) return;
	if (pointId == U.id) {
		// calculate angle
		double SUE = angleBetweenVect(vect(S, U), vect(E, U));
		double SDE = angleBetweenVect(vect(S, D), vect(E, D));
		if (SDE < 0 || SUE < 0) cout << "negative angle error, should not come\n";

		if ((SDE + SUE) <= pi) return;	// outside delaunay edge
		else {	// inside, so flip edge
			edgeList[edgeId].isValid = false;
			mapping[edgeList[edgeId].startP.id][edgeList[edgeId].endP.id] = -1;
			mapping[edgeList[edgeId].endP.id][edgeList[edgeId].startP.id] = -1;

			// flip edge create
			DelaunayEdge newEdge(U, D, S, E);
			newEdge.isValid = 1;
			newEdge.id = numberOfInputEdges++;
			edgeList.push_back(newEdge);
			mapping[newEdge.startP.id][newEdge.endP.id] = newEdge.id;
			mapping[newEdge.endP.id][newEdge.startP.id] = newEdge.id;

			// update old edge
			// SU->E , D
			if (edgeList[mapping[S.id][U.id]].uP.id == E.id) edgeList[mapping[S.id][U.id]].uP = D;
			else if (edgeList[mapping[S.id][U.id]].dP.id == E.id) {
				edgeList[mapping[S.id][U.id]].dP = D;
			}
			// SD->E, U
			if (edgeList[mapping[S.id][D.id]].uP.id == E.id) edgeList[mapping[S.id][D.id]].uP = U;
			else if (edgeList[mapping[S.id][D.id]].dP.id == E.id) {
				edgeList[mapping[S.id][D.id]].dP = U;
			}
			// ED->S, U
			if (edgeList[mapping[E.id][D.id]].uP.id == S.id) edgeList[mapping[E.id][D.id]].uP = U;
			else if (edgeList[mapping[E.id][D.id]].dP.id == S.id) {
				edgeList[mapping[E.id][D.id]].dP = U;
			}
			// EU->S, D
			if (edgeList[mapping[E.id][U.id]].uP.id == S.id) edgeList[mapping[E.id][U.id]].uP = D;
			else if (edgeList[mapping[E.id][U.id]].dP.id == S.id) {
				edgeList[mapping[E.id][U.id]].dP = D;
			}

			int edgeSD = mapping[S.id][D.id];
			int edgeED = mapping[E.id][D.id];
			//calling legalize
			legalizeEdge(pointId, edgeSD);
			legalizeEdge(pointId, edgeED);
		}
	}
	else if (D.id > 0 && pointId == D.id) {	// making sure D is not -1
		// calculate angle
		double SUE = angleBetweenVect(vect(S, U), vect(E, U));
		double SDE = angleBetweenVect(vect(S, D), vect(E, D));
		if (SDE < 0 || SUE < 0) cout << "negative angle error, should not come\n\n";

		if ((SDE + SUE) <= pi) return;	// outside delaunay edge
		else {	// inside, so flip edge
			edgeList[edgeId].isValid = false;
			mapping[edgeList[edgeId].startP.id][edgeList[edgeId].endP.id] = -1;
			mapping[edgeList[edgeId].endP.id][edgeList[edgeId].startP.id] = -1;

			// flip edge create
			DelaunayEdge newEdge(U, D, S, E);
			newEdge.isValid = 1;
			newEdge.id = numberOfInputEdges++;
			edgeList.push_back(newEdge);
			mapping[newEdge.startP.id][newEdge.endP.id] = newEdge.id;
			mapping[newEdge.endP.id][newEdge.startP.id] = newEdge.id;

			// update old edge
			// SU->E , D
			if (edgeList[mapping[S.id][U.id]].uP.id == E.id) edgeList[mapping[S.id][U.id]].uP = D;
			else if (edgeList[mapping[S.id][U.id]].dP.id == E.id) {
				edgeList[mapping[S.id][U.id]].dP = D;
			}
			// SD->E, U
			if (edgeList[mapping[S.id][D.id]].uP.id == E.id) edgeList[mapping[S.id][D.id]].uP = U;
			else if (edgeList[mapping[S.id][D.id]].dP.id == E.id) {
				edgeList[mapping[S.id][D.id]].dP = U;
			}
			// ED->S, U
			if (edgeList[mapping[E.id][D.id]].uP.id == S.id) edgeList[mapping[E.id][D.id]].uP = U;
			else if (edgeList[mapping[E.id][D.id]].dP.id == S.id) {
				edgeList[mapping[E.id][D.id]].dP = U;
			}
			// EU->S, D
			if (edgeList[mapping[E.id][U.id]].uP.id == S.id) edgeList[mapping[E.id][U.id]].uP = D;
			else if (edgeList[mapping[E.id][U.id]].dP.id == S.id) {
				edgeList[mapping[E.id][U.id]].dP = D;
			}

			int edgeSU = mapping[S.id][U.id];
			int edgeEU = mapping[E.id][U.id];
			//calling legalize
			legalizeEdge(pointId, edgeSU);
			legalizeEdge(pointId, edgeEU);
		}
	}

}
