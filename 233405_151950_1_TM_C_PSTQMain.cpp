// ..

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <conio.h>
#include <math.h>
#include<vector>
#include <list>
#include<windows.h>
//#include "RTree.h"
using namespace std;
typedef int ValueType;

class Node{  // Class to handle Node objects
public:
	int name; //Name index as used in the paper
	double x;  //Coordinate X
	double y;	//Coordinate Y
	void set(int n, double xx, double yy)
	{
		name = n;
		x = xx;
		y = yy;
	}
	void disp()
	{
		cout<<"n"<<name<<"   Location: ("<<x<<", "<<y<<")";
	}
};

class PDF{ // Class to handle Discrete PDFs (3 parts) as used in the paper
public:
	double v[3];	//PDF function values
	double p[3];	//PDF function probabilities
	void disp()
	{
		for (int i=0; i<3; i++)
		cout<< "("<<v[i]<<", "<<p[i]<<")"; 
	}
};

class EdgeTimes{ // Class to handle travel times in Edge
public:
	double t[3]; // Time Value
	double p[3]; // Probability
	void set(double distance, PDF pdf)
	{
		for (int i=0; i<3;i++)
		{
			p[i] = pdf.p[i]; // Probaabilities transferred 
			t[i] = distance / pdf.v[i];  //time=distance/speed
		}
	}
	void disp()
	{
		for (int i=0; i<3; i++)
		cout<< "("<<t[i]<<", "<<p[i]<<")"; 
	}
};

class PathTimes{ // Class to handle travel times in Path
public:
	double t[1024]; // Time Value
	double p[1024]; // Probability
};

class Edge{	// Class to handle Edge objects
public:
	int name; //Name index as used in the paper
	Node from;	//Starting Node
	Node to;	//End Node
	double distance;	//Edge Distance
	PDF pdf;
	EdgeTimes traveltime;
	bool walkingFromTo; //This flag is true when we walk from-->to (During traversal)
	// It is false when we walk from to-->from (During traversal)

	void changeDirection()
	{
		if (walkingFromTo)
			walkingFromTo = false;
		else
			walkingFromTo = true;
	}

	void set(int n, Node a, Node b)
	{
		name = n;
		from = a;
		to = b;
		distance = sqrt(double((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y)));	
		walkingFromTo = true;
	}
	void setpdf(double v1, double p1, double v2, double p2, double v3, double p3)
	{
		pdf.v[0] = v1; pdf.v[1] = v2; pdf.v[2] = v3;
		pdf.p[0] = p1; pdf.p[1] = p2; pdf.p[2] = p3;
		traveltime.set(distance, pdf);
	}
	void setpdfObject(PDF p)
	{
		pdf = p;
		traveltime.set(distance, pdf);
	}
	
	void disp()
	{
		if (-1 == name)
		{
			if ((-1 == from.name) && (-1 == to.name)) 
				cout<< "e (from,to)    Distance:"<<distance<<"\t\t"; 
			if ((-1 != from.name) && (-1 == to.name)) 
				cout<< "e ("<<from.name<<",to)    Distance:"<<distance<<"\t\t"; 
			if ((-1 == from.name) && (-1 != to.name)) 
				cout<< "e (from,"<<to.name<<")    Distance:"<<distance<<"\t\t"; 
			if ((-1 != from.name) && (-1 != to.name)) 
				cout<< "e ("<<from.name<<","<<to.name<<")    Distance:"<<distance<<"\t\t"; 
		}
		else
		{
			if ((-1 == from.name) && (-1 == to.name)) 
				cout<< "e"<<name<<"(from,to)    Distance:"<<distance<<"\t\t"; 
			if ((-1 != from.name) && (-1 == to.name)) 
				cout<< "e"<<name<<"("<<from.name<<",to)    Distance:"<<distance<<"\t\t"; 
			if ((-1 == from.name) && (-1 != to.name)) 
				cout<< "e"<<name<<"(from,"<<to.name<<")    Distance:"<<distance<<"\t\t"; 
			if ((-1 != from.name) && (-1 != to.name)) 
				cout<< "e"<<name<<"("<<from.name<<","<<to.name<<")    Distance:"<<distance<<"\t\t"; 
		}

		if (0 == (distance - int(distance)))
			cout<<"\t";
		pdf.disp();
		cout<<"\t";
		traveltime.disp();
	}
	void dispmin()
	{
		if (-1 == name)
		{
			if ((-1 == from.name) && (-1 == to.name)) 
				cout<< "e (from,to)"; 
			if ((-1 != from.name) && (-1 == to.name)) 
				cout<< "e ("<<from.name<<",to)"; 
			if ((-1 == from.name) && (-1 != to.name)) 
				cout<< "e (from,"<<to.name<<")"; 
			if ((-1 != from.name) && (-1 != to.name)) 
				cout<< "e ("<<from.name<<","<<to.name<<")"; 
		}
		else
		{
			if ((-1 == from.name) && (-1 == to.name)) 
				cout<< "e"<<name<<"(from,to)"; 
			if ((-1 != from.name) && (-1 == to.name)) 
				cout<< "e"<<name<<"("<<from.name<<",to)"; 
			if ((-1 == from.name) && (-1 != to.name)) 
				cout<< "e"<<name<<"(from,"<<to.name<<")"; 
			if ((-1 != from.name) && (-1 != to.name)) 
				cout<< "e"<<name<<"("<<from.name<<","<<to.name<<")"; 
		}
	}
	
};

class Facilities{ // Class to handle Facilities
public:
	int name; //Name index as used in the paper
	double x;  //Coordinate X
	double y;	//Coordinate Y
	int type;   // Type of Facility
	//  Type 1: Hotel 
	//  Type 2: Restaurant
	//  Type 3: Airport
	string typetext;
	Edge e;
	void set(int n, double xx, double yy, int t, Edge ee)
	{
		name = n;
		x = xx;
		y = yy;
		type = t;
		switch(t)
		{
		case 0: { typetext = "Query Pos."; break;} // Current Position
		case 1: { typetext = "Hotel     "; break;}
		case 2: { typetext = "Restaurant"; break;}
		case 3: { typetext = "Airport   "; break;}
		default: { typetext = "NoName    "; break;}
		};
		e.set(ee.name, ee.from, ee.to);
		e.setpdfObject(ee.pdf);
	}
	void disp()
	{
		cout<< "o"<<name<<"("<<typetext<<")  \t  Position: ("<<x<<", "<<y<<")"<<"\t";
		e.dispmin();
	}
};

class Path{ //Class to handle Path objects
public:
	int name; //Name index
	Facilities from;
	Facilities to;
	Edge ViaE[5]; // An array of Edges on the path 
	// Note that the Edge array ViaE does not include Edges 
	// where Facilities from and to are located
	// Maximum here is taken as 5, this can be later modified to 
	// Have a dynamic size
	int MaxNumViaE;
	int NumViaE; // Current meaningful size of ViaE array 

	Edge E[7]; // The complete Path
	// The complete path includes the following:
	// E[0]: Edge from Facility from to first ViaE
	// E[1..(NumViaE-1)]:  Via Edges
	// E[NumViaE+1] (or E[NumE-1]): Edge from ViaE[NumViaE-1] to Facility to
	int NumE;

	long numPossibleTimes;
	long numPossibleTimesCounter;
	PathTimes times;
	double Tmin, Tmax;
	double Pmin, Pmax;

	double PSTQ;
	void computeTimes()
	{
		numPossibleTimes = long (pow(float(3), float(NumE)));  
		numPossibleTimesCounter = 0;
		computeTimesRecursiveStep(1,0,1);

		sortTimesByProb();
		Pmin = times.p[numPossibleTimes-1];
		Pmax = times.p[0];

		sortTimesByValue();				
		Tmin = times.t[0];
		Tmax = times.t[numPossibleTimes-1];
	}

	void computeTimesRecursiveStep(int level, double preT, double preP)
	{
		if (level > NumE)
		{
			times.t[numPossibleTimesCounter] = preT;
			times.p[numPossibleTimesCounter] = preP;
			numPossibleTimesCounter++;
			
		}
		else
		{
			double newT1, newT2, newT3;
			double newP1, newP2, newP3;
			newT1 = double(preT + E[level-1].traveltime.t[0]);
			newT2 = double(preT + E[level-1].traveltime.t[1]);
			newT3 = double(preT + E[level-1].traveltime.t[2]);
			newP1 = double(preP * E[level-1].traveltime.p[0]);
			newP2 = double(preP * E[level-1].traveltime.p[1]);
			newP3 = double(preP * E[level-1].traveltime.p[2]);
			
			level++;
			computeTimesRecursiveStep(level, newT1, newP1);
			computeTimesRecursiveStep(level, newT2, newP2);
			computeTimesRecursiveStep(level, newT3, newP3);
		}
	}

	void sortTimesByValue()
	{
		for (int i=0; i< (numPossibleTimes-1);i++)
		{
			double min = times.t[i];
			int swapIndex = i;
			for (int j=(i+1); j<numPossibleTimes; j++)
			{
				if(times.t[j] < min)
				{
					min = times.t[j];
					swapIndex = j;
				}
			}
			double swapperT = times.t[i];
			double swapperP = times.p[i];
			times.t[i] = times.t[swapIndex];
			times.p[i] = times.p[swapIndex];
			times.t[swapIndex] = swapperT;
			times.p[swapIndex] = swapperP;
		}
	}

	void sortTimesByProb()
	{
		for (int i=0; i< (numPossibleTimes-1);i++)
		{
			double max = times.p[i];
			int swapIndex = i;
			for (int j=(i+1); j<numPossibleTimes; j++)
			{
				if(times.p[j] > max)
				{
					max = times.t[j];
					swapIndex = j;
				}
			}
			double swapperT = times.t[i];
			double swapperP = times.p[i];
			times.t[i] = times.t[swapIndex];
			times.p[i] = times.p[swapIndex];
			times.t[swapIndex] = swapperT;
			times.p[swapIndex] = swapperP;
		}
	}

	void displayTimes()
	{
		for(long i=0; i<numPossibleTimes; i++)
		{
			cout<<"\nTime Possibility #"<<i+1<<" : t "<<times.t[i]<<"\t p "<<times.p[i];
		}
	}

	void displayEdgeTimes()
	{
		for (int k=0; k< NumE; k++)
		{
			cout<<"\nEdge "<<k;
			for(long i=0; i<3; i++)
			{
				cout<<"\nTime Possibility #"<<i+1<<" : t "<<E[k].traveltime.t[i]<<"\t p "<<E[k].traveltime.p[i];
			}
		}
	}

	void displayPathEdges()
	{
		for(int i=0; i<NumE; i++)
		{
			cout<<"\n";
			E[i].dispmin();
		}
	}

	void displayPathNodes()
	{
		for(int i=0; i<NumE; i++)
		{
			if (-1 == E[i].from.name)
				cout<<"q";
			else
				cout<<"n"<<E[i].from.name;
			cout<<"-->";
		}
		cout<<"o";
	}

	void set(int n, Facilities f, Edge a, Edge b, Edge c, Edge d, Edge e, Facilities t, int num=5 )
	{
		name = n;
		from = f;
		ViaE[0] = a; ViaE[1] = b; ViaE[2] = c; ViaE[3] = d; ViaE[4] = e;
		to = t;
		NumViaE = num;
		MaxNumViaE = 5;
		NumE = num+2;
		
		Edge tempEdge;
		Node tempNode;
		tempNode.set(-1, f.x, f.y);
		if (ViaE[0].walkingFromTo)
		{
			tempEdge.set(-1, tempNode, ViaE[0].from);
			tempEdge.setpdfObject(f.e.pdf);
		}
		else
		{
			tempEdge.set(-1, tempNode, ViaE[0].to);
			tempEdge.setpdfObject(f.e.pdf);
		}
		
		E[0] = tempEdge;
		
		for (int i=0; i<NumViaE; i++)
			E[i+1] = ViaE[i];
		
		tempNode.set(-1, t.x, t.y);
		if (ViaE[NumViaE].walkingFromTo)
		{
			tempEdge.set(-1, ViaE[NumViaE-1].to, tempNode);
			tempEdge.setpdfObject(t.e.pdf);
		}
		else
		{
			tempEdge.set(-1, ViaE[NumViaE-1].from, tempNode);
			tempEdge.setpdfObject(t.e.pdf);
		}
		E[NumE-1] = tempEdge;
	}
	void set(int n, Facilities f, Facilities t)
	{
		name = n;
		from = f;
		to = t;
		NumViaE = 0;
		MaxNumViaE = 5;
		NumE = 1;

		Node tempNodef, tempNodet;
		tempNodef.set(-1, f.x, f.y);
		tempNodet.set(-1, t.x, t.y);
		Edge tempEdge;
		tempEdge.set(-1, tempNodef, tempNodet);
		tempEdge.setpdfObject(f.e.pdf);			//Dirty Line
		E[0]=tempEdge;		
	}
	bool addE(Edge e)
	{
		if (NumViaE <= MaxNumViaE)
		{
			if((1 == NumE) && (0 == NumViaE)) // Entering first ViaE (Double Check)
			{
				NumViaE++;
				ViaE[NumViaE-1] = e;

				NumE = NumE + 2; // The first entry of ViaE edge will give us 3 total edges

				Edge tempEdge;
				Node tempNode;
				tempNode.set(-1, from.x, from.y);
				if (ViaE[0].walkingFromTo)
				{
					tempEdge.set(-1, tempNode, ViaE[0].from);
					tempEdge.setpdfObject(from.e.pdf);
				}
				else
				{
					tempEdge.set(-1, tempNode, ViaE[0].to);
					tempEdge.setpdfObject(from.e.pdf);
				}
		
				E[0] = tempEdge;
		
				for (int i=0; i<NumViaE; i++)
					E[i+1] = ViaE[i];
		
				tempNode.set(-1, to.x, to.y);
				if (ViaE[NumViaE].walkingFromTo)
				{
					tempEdge.set(-1, ViaE[NumViaE-1].to, tempNode);
					tempEdge.setpdfObject(to.e.pdf);
				}
				else
				{
					tempEdge.set(-1, ViaE[NumViaE-1].from, tempNode);
					tempEdge.setpdfObject(to.e.pdf);
				}
				E[NumE-1] = tempEdge;
			}
			else
			{
				NumViaE++;
				NumE++;

				ViaE[NumViaE-1] = e;
				E[NumE-2] = e;

				Edge tempEdge;
				Node tempNode;

				tempNode.set(-1, to.x, to.y);
				if (ViaE[NumViaE].walkingFromTo)
				{
					tempEdge.set(-1, ViaE[NumViaE-1].to, tempNode);
					tempEdge.setpdfObject(to.e.pdf);
				}
				else
				{
					tempEdge.set(-1, ViaE[NumViaE-1].from, tempNode);
					tempEdge.setpdfObject(to.e.pdf);
				}
				
				E[NumE-1] = tempEdge;
			}

			return true;
		}
		else
			return false;
	}
};


double PSTQ(Path The, Path Other1, Path Other2)
{
	double Sigma = 0;
	for(long i=0; i<The.numPossibleTimes; i++)
	{
		double Iteration = The.times.p[i];

		double probSum1 = 0;
		for(int j=0; j< Other1.numPossibleTimes; j++)
			if (Other1.times.t[j] <= The.times.t[i] )
				probSum1 += Other1.times.p[j];
		Iteration *= (1-probSum1);

		double probSum2 = 0;
		for(int j=0; j< Other2.numPossibleTimes; j++)
			if (Other2.times.t[j] <= The.times.t[i] )
				probSum2 += Other2.times.p[j];
		Iteration *= (1-probSum2);
		
		Sigma += Iteration;
	}
	return Sigma;
}

double PSTQ_List(Path The, list<Path> Other)
{
	double Sigma = 0;
	for(long i=0; i<The.numPossibleTimes; i++)
	{
		double Iteration = The.times.p[i];

		list <Path> cOther = Other;
		for(int k=0; k<Other.size(); k++)
		{
			Path temp = cOther.back();
			cOther.pop_back();
			double probSum = 0;
			for(int j=0; j< temp.numPossibleTimes; j++)
			if (temp.times.t[j] <= The.times.t[i] )
				probSum += temp.times.p[j];
			Iteration *= (1-probSum);
		}

		Sigma += Iteration;
	}
	return Sigma;
}

struct Rect
{
  Rect()  {}

  Rect(int a_minX, int a_minY, int a_maxX, int a_maxY)
  {
    min[0] = a_minX;
    min[1] = a_minY;

    max[0] = a_maxX;
    max[1] = a_maxY;
  }


  int min[2];
  int max[2];
};

struct Rect rects[] =
{
  Rect(0, 0, 2, 2), // xmin, ymin, xmax, ymax (for 2 dimensional RTree)
  Rect(5, 5, 7, 7),
  Rect(8, 5, 9, 6),
  Rect(7, 1, 9, 2),
};

int nrects = sizeof(rects) / sizeof(rects[0]);

Rect search_rect(6, 4, 10, 6); // search will find above rects that this one overlaps


bool MySearchCallback(ValueType id, void* arg)
{
  cout << "Hit data rect " << id << "\n";
  return true; // keep going
}


void oldmain(bool getchflag = true)
{
	// Put this flag as false if you donot want to see the output in step by step manner.
	
	cout<<"PSTQ Proof-of-Concept Implementation\n\n";
	

	// Declaring Nodes
	Node *n;
	n = new Node[11];  // 11 Nodes as defined in the paper
	int NumNode = 11;
	// Defining Nodes
	n[0].set(1,5,1);	n[1].set(2,3,2);	n[2].set(3,1,3);	n[3].set(4,2,5);	n[4].set(5,4,5);
	n[5].set(6,7,5);	n[6].set(7,9,5);	n[7].set(8,9,3);	n[8].set(9,9,0);	n[9].set(10,8,0);
	n[10].set(11,7,2);

	// Declaring Edges
	Edge *e;
	e = new Edge[14]; //14 Edges as defined in the paper
	int NumEdge = 14;
	// Defining Edges and PDFs for the Edges
	e[0].set(1,n[0],n[1]);		e[0].setpdf(30,0.5,35,0.3,36,0.2);
	e[1].set(2,n[1],n[2]);		e[1].setpdf(45,0.6,48,0.2,50,0.2);
	e[2].set(3,n[2],n[3]);		e[2].setpdf(50,0.4,52,0.3,55,0.3);
	e[3].set(4,n[3],n[4]);		e[3].setpdf(45,0.3,46,0.3,50,0.4);
	e[4].set(5,n[1],n[4]);		e[4].setpdf(55,0.6,58,0.3,60,0.1);
	e[5].set(6,n[4],n[5]);		e[5].setpdf(47,0.2,48,0.3,50,0.5);
	e[6].set(7,n[5],n[6]);		e[6].setpdf(48,0.3,49,0.3,50,0.4);
	e[7].set(8,n[5],n[10]);		e[7].setpdf(30,0.7,32,0.2,35,0.1);
	e[8].set(9,n[6],n[7]);		e[8].setpdf(50,0.3,51,0.3,52,0.4);
	e[9].set(10,n[0],n[10]);	e[9].setpdf(36,0.1,38,0.4,40,0.5);
	e[10].set(11,n[7],n[8]);	e[10].setpdf(40,0.6,42,0.3,45,0.1);
	e[11].set(12,n[8],n[9]);	e[11].setpdf(40,0.8,42,0.1,45,0.1);
	e[12].set(13,n[9],n[10]);	e[12].setpdf(55,0.3,58,0.3,60,0.4);
	e[13].set(14,n[7],n[10]);	e[13].setpdf(40,0.4,42,0.3,45,0.3);

	// Declaring Facilities
	Facilities *o;
	o = new Facilities[5]; //5 Facilities as defined in the paper
	int NumFacilitie = 5;
	// Defining Facilities
	o[0].set(1,6,5,1,e[5]);
	o[1].set(2,7,4,2,e[7]);
	o[2].set(3,7.5,1,3,e[12]);
	o[3].set(4,2,2.5,1,e[1]);
	o[4].set(5,9,1,1,e[10]);


	// Defining The Query Position (In Paper this is Position q)
	// We use the Query position as a "Facility" object as that too
	// is a location, and exists on a edge
	Facilities q;
	q.set(0, 6, 1.5, 0,e[9]);
	// Query position q (Thus type=0) with location(6, 1.5) and on Edge e(1,11)

	// Declaring Paths form q to the three Hotels (as defined in the paper)
	Path *path;
	path = new Path[3];
	int NumPaths = 3;

	// Path to o1 Hotel
	path[0].set(1, q, o[0]);
	Edge tmp1;	tmp1.set(e[7].name, e[7].to, e[7].from);	tmp1.setpdfObject(e[7].pdf);
	// As e[7] is from 6 to 11 and we are walking from 11 to 6
	path[0].addE(tmp1);

	// Path to o4 Hotel
	path[1].set(2, q, o[3]);
	path[1].addE(e[0]);

	// Path to o5 Hotel
	path[2].set(3, q, o[4]);
	Edge tmp2;	tmp2.set(e[13].name, e[13].to, e[13].from);	tmp2.setpdfObject(e[13].pdf);
	// As e[13] is from 8 to 11 and we are walking from 11 to 8
	path[2].addE(tmp2);


	//void set(int n, Facilities f, Edge a, Edge b, Edge c, Edge d, Edge e, Facilities t, int num=5 )


	cout<<"\n------------------------Nodes:";
	for (int i=0; i< NumNode; i++)
	{
		cout<<"\n";
		n[i].disp();

	}
	if (getchflag) getch();
	cout<<"\n------------------------Edges:              PDFs                         Times";
	for (int i=0; i< NumEdge; i++)
	{
		cout<<"\n";
		e[i].disp();
	}
	if (getchflag) getch();
	cout<<"\n------------------------Facilities:";
	for (int i=0; i< NumFacilitie; i++)
	{
		cout<<"\n";
		o[i].disp();
	}
	if (getchflag) getch();
	cout<<"\n------------------------Paths from q to Hotels:";
	for (int i=0; i< NumPaths; i++)
	{
		cout<<"\n--------------------------- Path "<< path[i].name;
		cout<<"[";
		path[i].displayPathNodes();
		cout<<"]";
		//path[i].displayPathEdges();
		cout<<"\n------------Times:\n";
		path[i].computeTimes();
		path[i].displayEdgeTimes();
		cout<<"\nFull Path";
		path[i].displayTimes();
		cout<<"\nT+ : "<< path[i].Tmin <<"  T- : "<< path[i].Tmax;
		cout<<"\tP- : "<< path[i].Pmin <<"  P+ : "<< path[i].Pmax;
	}
	cout<<"\n--NOTE THAT THESE VALUES ARE A BIT DIFFERENT THAN THE ONES SHOWN IN THE PAPER AS IN PAPER THE PROBABILITIES WHICH ARE VERY SIMILAR ARE JOINED TOGETHER AS ONE VALUE";
	/*   //A Test Path 
	Path p = path[0];
	p.addE(tmp2);
	p.addE(tmp2);
	p.displayPathNodes();
	p.computeTimes();
	p.displayTimes();
	*/
	if (getchflag) getch();
	//Computing PSTQ for the three paths
	path[0].PSTQ = PSTQ(path[0], path[1], path[2]);
	path[1].PSTQ = PSTQ(path[1], path[0], path[2]);
	path[2].PSTQ = PSTQ(path[2], path[1], path[0]);
	cout<<"\n------------------------PSTQ for the three Paths:";
	double TotalPSTQ = 0;
	for(int i=0; i<3; i++)
	{
		cout<<"\nPath "<<i+1<<"\tPSTQ:"<<path[i].PSTQ;
		TotalPSTQ += path[i].PSTQ;
	}
	cout<<"\nTotal  \tPSTQ:"<<TotalPSTQ;
	cout<<"\n--NOTE THAT THESE VALUES ARE A BIT DIFFERENT THAN THE ONES SHOWN IN THE PAPER AS IN PAPER THE PROBABILITIES WHICH ARE VERY SIMILAR ARE JOINED TOGETHER AS ONE VALUE";

	double alpha = 0.80;   //80%
	cout<<"\nOutput for Query PSTQ(alpha = 80%)";
	for(int i=0; i<3; i++)
	{
		cout<<"\nStep Count:"<<i+1;
		if (path[i].PSTQ >= alpha )
		cout<<"\nPath "<<i+1<<"\tPSTQ:"<<path[i].PSTQ;		
	}

	if (getchflag) getch();
	//Time Bound Pruning
	cout<<"\nOutput for Query PSTQ(alpha = 80%) [Time Bound Pruning]";
	double MIN = 0;
	for(int i=0; i<3; i++)
	{
		if (0 == MIN) MIN = path[0].Tmin;;
		cout<<"\nStep Count [Time Bound Pruning]:"<<i+1;
		if (path[i].Tmin <= MIN)
		{
			MIN = path[i].Tmin;
			if (path[i].PSTQ >= alpha )
			cout<<"\nPath "<<i+1<<"\tPSTQ:"<<path[i].PSTQ;
		}
		else
			cout<<"\tSkip Step";
	}
	if (getchflag) getch();
	//Probabilistic Bound Pruning
	
	double beta = (1-alpha);   //	Assuming the conditions for beta holds  // You may put your own value here for beta

	cout<<"\nOutput for Query PSTQ(Alpha = 80%, Beta = "<<beta<<") [Time Bound Pruning]";
	double MAX = 0;
	for(int i=0; i<3; i++)
	{
		if (0 == MAX) MAX = path[0].Tmin;
		cout<<"\nStep Count [Probabilistic Bound Pruning]:"<<i+1;
		if ( path[i].Tmin > (MAX* beta))
		{
			MAX = path[i].Tmax;
			if (path[i].PSTQ >= alpha )
			cout<<"\nPath "<<i+1<<"\tPSTQ:"<<path[i].PSTQ;
		}
		else
			cout<<"\tSkip Step";
	}
	if (getchflag) getch();
	cout<<"\n\n";
	cout<<"\n//////////////////////////R* Tree Optimization on Pre Computed Data//////////////////////////////////";
	

	/*
	Here is the RTree class and example (Modified version of the R Tree standard header available on GITHUB)
	Try and understand this. You should be able to fill in the Nodal and Leaf values to the example to implement it
	for PSTQ.
	*/
	/*
	
	typedef RTree<ValueType, int, 2, float> MyTree;
  MyTree tree;

  int i, nhits;
  cout << "nrects = " << nrects << "\n";

  for(i=0; i<nrects; i++)
  {
    tree.Insert(rects[i].min, rects[i].max, i); // Note, all values including zero are fine in this version
  }

  nhits = tree.Search(search_rect.min, search_rect.max, MySearchCallback, NULL);

  cout << "Search resulted in " << nhits << " hits\n";

  // Iterator test
  int itIndex = 0;
  MyTree::Iterator it;
  for( tree.GetFirst(it);
       !tree.IsNull(it);
       tree.GetNext(it) )
  {
    int value = tree.GetAt(it);

    int boundsMin[2] = {0,0};
    int boundsMax[2] = {0,0};
    it.GetBounds(boundsMin, boundsMax);
    cout << "it[" << itIndex++ << "] " << value << " = (" << boundsMin[0] << "," << boundsMin[1] << "," << boundsMax[0] << "," << boundsMax[1] << ")\n";
  }

  // Iterator test, alternate syntax
  itIndex = 0;
  tree.GetFirst(it);
  while( !it.IsNull() )
  {
    int value = *it;
    ++it;
    cout << "it[" << itIndex++ << "] " << value << "\n";
  }
  */
  getch();

   
}

void GRAPHHER( bool flag);

int main()
{
	
	bool flag = true;
	while(flag)
	{
		char a,b;
		cout<<"\n***Start Menu***********************************";
		cout<<"\n[1] View Paper implementation with values given in the paper.";
		cout<<"\n[2] Make a new query on the map.";
		cout<<"\n[c] Clear Screen.";
		cout<<"\n[0] Exit";
		cout<<"\n***********************************************     :\>";
		
		cin>>a;
		switch(a)
		{
		case '1': { oldmain(true); break;}
		case '2': { 
			bool flag2 = true;
			while(flag2)
			{
			cout<<"\n***Map Menu***********************************";
			cout<<"\n[1] Take the map given in the paper.";
			cout<<"\n[2] Start With a blank map.";
			cout<<"\n[c] Clear Screen.";
			cout<<"\n[0] Exit";
			cout<<"\n***********************************************     :\>";
			cin>>b;
			switch(b)
			{
			case '1':{ GRAPHHER(true);break;}
			case '2':{ GRAPHHER(false); break;}
			case 'c':
			case 'C': { cout<<string(100,'\n'); break;}
			case '0': {flag2 = false; break;}
			default: break;
			}

			}

			break;};
		case 'c':
		case 'C': { cout<<string(100,'\n'); break;}
		case '0': {flag = false; break;}
		default: break;
		}
	}
	
	cin.get();
	return 0;
}

void drawLineDDA(int xa, int ya, int xb, int yb, COLORREF COLOR, HDC mydc){

   int dx = xb - xa, dy = yb - ya, steps, k;
   float xIncrement, yIncrement, x = xa, y = ya;
   if(abs(dx) > abs(dy)) steps = abs(dx);
   else steps = abs(dy);
   xIncrement = dx / (float) steps;
   yIncrement = dy / (float) steps;
   
   for(int k = 0; k < steps; k++){
	   x += xIncrement;
	   y += yIncrement;
	   SetPixel(mydc,x,y,COLOR);
   }
}


void boxaround(int x, int y, int s, COLORREF COLOR, HDC mydc)
{
	int t = s/2;
	drawLineDDA(x-t,y-t,x+t,y-t, COLOR, mydc);
	drawLineDDA(x+t,y-t,x+t,y+t, COLOR, mydc);
	drawLineDDA(x+t,y+t,x-t,y+t, COLOR, mydc);
	drawLineDDA(x-t,y+t,x-t,y-t, COLOR, mydc);
}

void deltaaround(int x, int y, int s, COLORREF COLOR, HDC mydc)
{
	int t = s/2;
	drawLineDDA(x,y-t,x+t,y+t, COLOR, mydc);
	drawLineDDA(x+t,y+t,x-t,y+t, COLOR, mydc);
	drawLineDDA(x-t,y+t,x,y-t, COLOR, mydc);
}

/*
void lableprint( int x, int y, char c, HWND myconsole, int sx = 79, int sy = 49) {

    // Set up the handles for reading/writing:
    //HANDLE wHnd = GetStdHandle(STD_OUTPUT_HANDLE);
    //HANDLE rHnd = GetStdHandle(STD_INPUT_HANDLE);

    // Change the window title:
   // SetConsoleTitle(TEXT("Win32 Console Control Demo"));

    // Set up the required window size:
    SMALL_RECT windowSize = {0, 0, sx, sy};
    
    // Change the console window size:
   // SetConsoleWindowInfo(wHnd, TRUE, &windowSize);
    
    // Create a COORD to hold the buffer size:
    //COORD bufferSize = {80, 50};

    // Change the internal buffer size:
    //SetConsoleScreenBufferSize(wHnd, bufferSize);

    // Set up the character:
    CHAR_INFO letter;
    letter.Char.AsciiChar = c;
    
    letter.Attributes = 
        FOREGROUND_RED | FOREGROUND_INTENSITY |
        BACKGROUND_RED | BACKGROUND_GREEN | BACKGROUND_INTENSITY;

    // Set up the positions:
    COORD charBufSize = {1,1};
    COORD characterPos = {x,y};
    SMALL_RECT writeArea = {0,0,sx,sy}; 

    // Write the character:
    WriteConsoleOutputA(myconsole, &letter, charBufSize, characterPos, &writeArea);

    // Move the cursor down a row so we can see the character:
    printf("\n");

}
*/


void PrintGraphher(list<Node> n, list<Edge> e, list<Facilities> o, Facilities q, list<Path> p, Path P, int delay = 1000)
{
	system("cls");
	//Get a console handle
    HWND myconsole = GetConsoleWindow();
	RECT r;
	GetWindowRect(myconsole, &r); //stores the console's current dimensions
	MoveWindow(myconsole, r.left, r.top, 1500, 800, TRUE); // 800 width, 100 height
    //Get a handle to device context
    HDC mydc = GetDC(myconsole);
	

	//int DrawText(myconsole,"A",1,);





	//Graphic Origin
	int XO = 500;
	int YO = 500;
	int S = 50;  //ScaleUP
	int SY = -50;

	//Choose any color
    COLORREF COLOR= RGB(255,255,255); 
	COLORREF BLUE= RGB(0,0,255); 
	COLORREF RED= RGB(255,0,0); 
	COLORREF GREEN= RGB(0,255,0); 
	COLORREF YELLOW= RGB(255,255,0); 
	COLORREF AXIS= RGB(100,100,0); 
	COLORREF GRID= RGB(0,50,50); 

	// Drawing Axis
	drawLineDDA(XO,YO,XO + (S * 12),YO, AXIS, mydc);
	drawLineDDA(XO,YO,XO,YO + (SY * 8), AXIS, mydc);

	//Drawing Grid
	for(int i=1; i<=8; i++)
	drawLineDDA(XO,YO+ (SY * i),XO + (S * 12),YO+ (SY * i), GRID, mydc);
	for(int i=1; i<=12; i++)
	drawLineDDA(XO+ (S * i),YO,XO + (S * i),YO+ (SY * 8), GRID, mydc);



	



	while(n.empty() != true)
	{
		Node temp = n.front();
		n.pop_front();
		SetPixel(mydc,XO + (S * temp.x),YO + (SY * temp.y),COLOR);
		boxaround(XO + (S * temp.x),YO + (SY * temp.y), S/4, BLUE, mydc);	
	}

	while(o.empty() != true)
	{
		Node temp;
		temp.set(o.front().name, o.front().x, o.front().y);
		o.pop_front();
		SetPixel(mydc,XO + (S * temp.x),YO + (SY * temp.y),COLOR);
		deltaaround(XO + (S * temp.x),YO + (SY * temp.y), S/4, YELLOW, mydc);	
	}

	Node temp;
	temp.set(q.name, q.x, q.y);
	SetPixel(mydc,XO + (S * temp.x),YO + (SY * temp.y),GREEN);
	deltaaround(XO + (S * temp.x),YO + (SY * temp.y), S/4, RED, mydc);
	deltaaround(XO + (S * temp.x),YO + (SY * temp.y), S/2, RED, mydc);
	deltaaround(XO + (S * temp.x),YO + (SY * temp.y), S/1.5, RED, mydc);


	while(e.empty() != true)
	{
		Edge temp = e.front();
		e.pop_front();
		drawLineDDA(XO + (S * temp.from.x),YO + (SY * temp.from.y),XO + (S * temp.to.x),YO + (SY * temp.to.y), COLOR, mydc);		
	}


	while(p.empty() != true)
	{
		Path temp = p.front();
		p.pop_front();

		for(int j=0; j<temp.NumE; j++)
		{
			Edge tempe = temp.E[j];
			drawLineDDA(XO + (S * tempe.from.x),YO + (SY * tempe.from.y),XO + (S * tempe.to.x),YO + (SY * tempe.to.y), YELLOW, mydc);	
			Beep( 750, 300 );
			Sleep(delay);
		}
		for(int j=0; j<temp.NumE; j++)
		{
			Edge tempe = temp.E[j];
			drawLineDDA(XO + (S * tempe.from.x),YO + (SY * tempe.from.y),XO + (S * tempe.to.x),YO + (SY * tempe.to.y), COLOR, mydc);	
		}
		
	}

	//Winning Path
	for(int j=0; j<P.NumE; j++)
	{
		Edge tempe = P.E[j];
		drawLineDDA(XO + (S * tempe.from.x),YO + (SY * tempe.from.y),XO + (S * tempe.to.x),YO + (SY * tempe.to.y), RED, mydc);	
		drawLineDDA(XO + (S * tempe.from.x),YO + (SY * tempe.from.y)+1,XO + (S * tempe.to.x),YO + (SY * tempe.to.y)+1, RED, mydc);	
		drawLineDDA(XO + (S * tempe.from.x),YO + (SY * tempe.from.y)-1,XO + (S * tempe.to.x),YO + (SY * tempe.to.y)-1, RED, mydc);	
		Beep( 1250, 300 );
		Sleep(delay/2);
	}
	
	/*
    //Draw pixels
	
	for(int i=500; i<1500; i++)
		for(int j=0; j<1500; j++)
			SetPixel(mydc,i,j,COLOR);

	drawLineDDA(10,10,350,350, BLUE, mydc);
	*/
	//lableprint( 150, 150, 'A', myconsole, 79,49);
	
	ReleaseDC(myconsole, mydc);
    //system("pause");
	
}



void GRAPHHER( bool flag)
{
	
	
	list<Node> n;
	Node tempnode;
	Node *nodes;
	list<Edge> e;
	Edge tempedge;
	Edge *edges;
	list<Facilities> o;
	Facilities tempfacility;
	Facilities *facilities;
	Facilities q;	
	list<Path> p;
	Path temppath;
	Path *paths;
	Path P;

	if(flag){
	tempnode.set(1,5,1);	n.push_back(tempnode);
	tempnode.set(2,3,2);	n.push_back(tempnode);
	tempnode.set(3,1,3);	n.push_back(tempnode);
	tempnode.set(4,2,5);	n.push_back(tempnode);
	tempnode.set(5,4,5);	n.push_back(tempnode);
	tempnode.set(6,7,5);	n.push_back(tempnode);
	tempnode.set(7,9,5);	n.push_back(tempnode);
	tempnode.set(8,9,3);	n.push_back(tempnode);
	tempnode.set(9,9,0);	n.push_back(tempnode);
	tempnode.set(10,8,0);	n.push_back(tempnode);
	tempnode.set(11,7,2);	n.push_back(tempnode);
	
	
	nodes = new Node[n.size()];
	list<Node> tt = n;
	for(int i=0; i< n.size(); i++)
	{
		nodes[i] = tt.front();
		tt.pop_front();
	}
	
	
	
	// Defining Edges and PDFs for the Edges
	tempedge.set(1,nodes[0],nodes[1]);		tempedge.setpdf(30,0.5,35,0.3,36,0.2);		e.push_back(tempedge);
	tempedge.set(2,nodes[1],nodes[2]);		tempedge.setpdf(45,0.6,48,0.2,50,0.2);		e.push_back(tempedge);
	tempedge.set(3,nodes[2],nodes[3]);		tempedge.setpdf(50,0.4,52,0.3,55,0.3);		e.push_back(tempedge);
	tempedge.set(4,nodes[3],nodes[4]);		tempedge.setpdf(45,0.3,46,0.3,50,0.4);		e.push_back(tempedge);
	tempedge.set(5,nodes[1],nodes[4]);		tempedge.setpdf(55,0.6,58,0.3,60,0.1);		e.push_back(tempedge);
	tempedge.set(6,nodes[4],nodes[5]);		tempedge.setpdf(47,0.2,48,0.3,50,0.5);		e.push_back(tempedge);
	tempedge.set(7,nodes[5],nodes[6]);		tempedge.setpdf(48,0.3,49,0.3,50,0.4);		e.push_back(tempedge);
	tempedge.set(8,nodes[5],nodes[10]);		tempedge.setpdf(30,0.7,32,0.2,35,0.1);		e.push_back(tempedge);
	tempedge.set(9,nodes[6],nodes[7]);		tempedge.setpdf(50,0.3,51,0.3,52,0.4);		e.push_back(tempedge);
	tempedge.set(10,nodes[0],nodes[10]);	tempedge.setpdf(36,0.1,38,0.4,40,0.5);		e.push_back(tempedge);
	tempedge.set(11,nodes[7],nodes[8]);	    tempedge.setpdf(40,0.6,42,0.3,45,0.1);		e.push_back(tempedge);
	tempedge.set(12,nodes[8],nodes[9]);	    tempedge.setpdf(40,0.8,42,0.1,45,0.1);		e.push_back(tempedge);
	tempedge.set(13,nodes[9],nodes[10]);	tempedge.setpdf(55,0.3,58,0.3,60,0.4);		e.push_back(tempedge);
	tempedge.set(14,nodes[7],nodes[10]);	tempedge.setpdf(40,0.4,42,0.3,45,0.3);		e.push_back(tempedge);

	
	edges = new Edge[e.size()];
	list<Edge> ttt = e;
	for(int i=0; i< e.size(); i++)
	{
		edges[i] = ttt.front();
		ttt.pop_front();
	}
	
	

	// Defining Facilities
	tempfacility.set(1,6,5,1,edges[5]);			o.push_back(tempfacility);
	tempfacility.set(2,7,4,2,edges[7]);			o.push_back(tempfacility);
	tempfacility.set(3,7.5,1,3,edges[12]);		o.push_back(tempfacility);
	tempfacility.set(4,2,2.5,1,edges[1]);		o.push_back(tempfacility);
	tempfacility.set(5,9,1,1,edges[10]);		o.push_back(tempfacility);
	
	
	facilities = new Facilities[o.size()];
	list<Facilities> tttt = o;
	for(int i=0; i< o.size(); i++)
	{
		facilities[i] = tttt.front();
		tttt.pop_front();
	}

	q.set(0, 6, 1.5, 0,edges[9]);


	// Declaring Paths form q to the three Hotels (as defined in the paper)
	
	
	// Path to o1 Hotel
	temppath.set(1, q, facilities[0]);
	Edge tmp1;	tmp1.set(edges[7].name, edges[7].to, edges[7].from);	tmp1.setpdfObject(edges[7].pdf);
	// As e[7] is from 6 to 11 and we are walking from 11 to 6
	temppath.addE(tmp1);
	p.push_back(temppath);

	// Path to o4 Hotel
	temppath.set(2, q, facilities[3]);
	temppath.addE(edges[0]);
	p.push_back(temppath);

	// Path to o5 Hotel
	temppath.set(3, q, facilities[4]);
	Edge tmp2;	tmp2.set(edges[13].name, edges[13].to, edges[13].from);	tmp2.setpdfObject(edges[13].pdf);
	// As e[7] is from 6 to 11 and we are walking from 11 to 6
	temppath.addE(tmp2);
	p.push_back(temppath);

	
	paths = new Path[3];
	list<Path> ttttt = p;
	for(int i=0; i< 3; i++)
	{
		paths[i] = ttttt.front();
		ttttt.pop_front();
	}
	P = paths[1];  //Winning Path
	}
	


	bool flag3 = true;

	while(flag3)
	{
		char c;
		int querytarget = 1;
		PrintGraphher(n,e,o,q,p,P, 100);
		for(int i=0; i< n.size(); i++)
			cout<<"\n Node ID: "<<i+1<<" Location("<<nodes[i].x<<", "<<nodes[i].y<<")";
		cout<<"\n";
		for(int i=0; i< e.size(); i++)
			cout<<"\n Edge ID: "<<i+1<<" Node "<<edges[i].from.name + 1<<" -> "<<edges[i].to.name + 1<<")";

		cout<<"\n***Work Menu***************";
		cout<<"\n[1] Add Node";
		cout<<"\n[2] Add Edge";
		cout<<"\n[3] Add Facility";
		cout<<"\n[4] Add New Query Point";
		cout<<"\n[c] Clear Screen.";
		cout<<"\n[0] Exit";
		cout<<"\n*************************** :\>";
		cin>>c;
		switch(c)
		{
		case '1':{ 
				Node temp;
				double x, y;
				cout<<"\nNode.x: ";
				cin>>x;
				cout<<"\nNode.y: ";
				cin>>y;
				temp.set(n.size(),x,y);
				n.push_back(temp);				
				nodes = new Node[n.size()];
				list<Node> tt = n;
				for(int i=0; i< n.size(); i++)
				{
					nodes[i] = tt.front();
					tt.pop_front();
				}
				
			break;}
		case '2':{
				Edge temp;
				int x, y;
				double v1,v2,v3,p1,p2,p3;
				cout<<"\nEdge.From Node: ";
				cin>>x;
				cout<<"\nEdge.To Node: ";
				cin>>y;
				cout<<"\nSpeed 1: "; cin>>v1;	cout<<"\nProbability 1: ";	cin>>p1;
				cout<<"\nSpeed 2: "; cin>>v2;	cout<<"\nProbability 2: ";	cin>>p2;
				cout<<"\nSpeed 3: "; cin>>v3;	cout<<"\nProbability 3: ";	cin>>p3;

				temp.set(e.size(),nodes[x-1],nodes[y-1]);
				temp.setpdf(v1,p1,v2,p2,v3,p3);
				e.push_back(temp);			
				edges = new Edge[e.size()];
				list<Edge> ttt = e;
				for(int i=0; i< e.size(); i++)
				{
					edges[i] = ttt.front();
					ttt.pop_front();
				}
				edges[e.size()-1].disp();
				Beep( 1250, 300 );	Sleep(5000);
			break;}
		case '3':{ 
				Facilities temp;
				double x, y;	int t; int choice;
				cout<<"\nFacility.x: ";
				cin>>x;
				cout<<"\nFacility.y: ";
				cin>>y;
				cout<<"\nFacility on Edge ID: ";
				cin>>choice;
				cout<<"\nFacility.Type [1 Hotel 2 Restaurant 3: Airport]: ";
				cin>>t;
				temp.set(o.size(),x,y,t,edges[choice]);
				o.push_back(temp);
				facilities = new Facilities[o.size()];
				list<Facilities> tttt = o;
				for(int i=0; i< o.size(); i++)
				{
					facilities[i] = tttt.front();
					tttt.pop_front();
				}
			
			break;}
		case '4':{ 
				Facilities temp;
				double x, y;	int t; int choice;
				cout<<"\nQuery Point.x: ";
				cin>>x;
				cout<<"\nQuery Point.y: ";
				cin>>y;
				cout<<"\nQuery Point on Edge ID: ";
				cin>>choice;			
				q.set(0, x, y, 0,edges[9]);
				while(p.empty() != true)
					p.pop_back();
				paths = NULL;
				cout<<"\nQuery Target: [1 Hotel 2 Restaurant 3: Airport]";
				cin>>querytarget;

				for(int i=0; i<o.size(); i++)
					if(facilities[i].type == querytarget)
						{
							cout<<"\n Candidate Facilities: ID "<<facilities[i].name +1;
							facilities[i].disp();
						}
				
				for(int i=0; i<o.size(); i++)
					if(facilities[i].type == querytarget)
					{
						cout<<"\n Enter Path for Candidate "<<i+1<<"  ";
						Path temp;
						temp.set(p.size(),q, facilities[i]);
						cout<<"\n";
						temp.displayPathNodes();
						int input = -1;
						while(input != 0)
						{
							cout<<"\nAdd intermediate Edge ID [0 for finish]:";
							cin>>input;
							if(input != 0)
								temp.addE(edges[input-1]);
							temp.displayPathEdges();
						}
						temp.computeTimes();
						temp.sortTimesByProb();
						p.push_back(temp);						
					}
				
				paths = new Path[p.size()];
				list<Path> ttttt = p;
				for(int i=0; i< p.size(); i++)
				{
					paths[i] = ttttt.front();
					ttttt.pop_front();
				}

				for(int i=0; i< p.size(); i++)
				{
						list<Path> otherpaths;
						for(int m=0; m< p.size(); m++)
						{
							if(m != i)
							{
								otherpaths.push_back(paths[m]);
							}
						}
					paths[i].PSTQ = PSTQ_List(paths[i], otherpaths);
					cout<<"\n";
					
				}
				double TotalPSTQ = 0;
				int winPath = 0;
				double winPSTQ = paths[0].PSTQ;
				for(int i=0; i<p.size(); i++)
				{
					cout<<"\nPath "<<i+1<<"\tPSTQ:"<<paths[i].PSTQ;
					TotalPSTQ += paths[i].PSTQ;

					if( winPSTQ < paths[i].PSTQ)
					{
						winPath = i;
						winPSTQ = paths[i].PSTQ;
					}

				}
				cout<<"\nTotal  \tPSTQ:"<<TotalPSTQ;
				cout<<"\n Winning path ID: "<< winPath+1;
				P = paths[winPath];
				Beep( 1250, 300 );	Sleep(5000);
				system("Pause");
				cout<<string(100,'\n');
			break;}
		case 'c':
		case 'C': { cout<<string(100,'\n'); break;}
		case '0': {flag3 = false; break;}
		default: break;
		}
	}


}




void xmain() {

    // Set up the handles for reading/writing:
    HANDLE wHnd = GetStdHandle(STD_OUTPUT_HANDLE);
    HANDLE rHnd = GetStdHandle(STD_INPUT_HANDLE);

    // Change the window title:
    SetConsoleTitle(TEXT("Win32 Console Control Demo"));

    // Set up the required window size:
    SMALL_RECT windowSize = {0, 0, 79, 49};
    
    // Change the console window size:
    SetConsoleWindowInfo(wHnd, TRUE, &windowSize);
    
    // Create a COORD to hold the buffer size:
    COORD bufferSize = {80, 50};

    // Change the internal buffer size:
    SetConsoleScreenBufferSize(wHnd, bufferSize);

    // Set up the character buffer:
    CHAR_INFO consoleBuffer[80 * 50];

    // We'll fill the console buffer with random data:
    //for (int y = 0; y < 50; ++y) {
      //  for (int x = 0; x < 80; ++x) {
        
            // An ANSI character is in the range 0-255,
            // so use % to keep it in this range.
            consoleBuffer[10 + 80 * 20].Char.AsciiChar = 'A';

            // The colour is also in the range 0-255,
            // as it is a pair of 4-bit colours.
            //consoleBuffer[10 + 80 * 20].Attributes = rand() % 256;
        //}
    //}

    // Set up the positions:
    COORD charBufSize = {80,50};
    COORD characterPos = {0,0};
    SMALL_RECT writeArea = {0,0,79,49}; 

    // Write the characters:
    WriteConsoleOutputA(wHnd, consoleBuffer, charBufSize, characterPos, &writeArea);

    // Move the cursor down a row so we can see the character:
    printf("\n");
	cin.get();
}