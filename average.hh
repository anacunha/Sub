// average.hh
//
// Class that performs the averaging step of a split/average subdivision
// sequence.

#ifndef AVERAGE_H
#define AVERAGE_H

class Cell;
class Vertex;
class Edge;

class Average
{
public:
	virtual void operator()( Cell *cell ) = 0;	// Perform averaging step on Cell
	virtual void applyEvaluation(Cell *cell) = 0;
};

class AvgNOOP : public Average
{
public:
	void operator()( Cell * );
	void applyEvaluation(Cell *cell);
	static void genNormal( Vertex *v );
	static void genNormals( Cell *c );
};

class AvgAdHoc : public Average
{
	bool interpolating;
public:
	AvgAdHoc( bool i = false )
		: interpolating( i )
	{
	}

	void operator()( Cell *cell );
	void applyEvaluation(Cell *cell);
	static void average( Vertex *v );
};

class AvgEval : public Average
{
	bool interpolating; // loop's OR butterfly
public:
	AvgEval(bool i = false)
		: interpolating(i)
	{
	}

	void operator()(Cell *cell);
	void applyEvaluation(Cell *cell);
	static void genNormal(Vertex *v);
	static void genNormals(Cell *cell);
	static void average(Vertex *v);
	static void butterfly(Vertex *v);
	static void extraordinary(Vertex *v, Vertex *v1, Edge *start, Edge *e, int n);
	static void evaluate(Vertex *v);
	static double tau(double i, double n);
};

#endif
