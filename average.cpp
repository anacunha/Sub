// average.cc - class that performs average step of subdivision

#include <assert.h>
#include <math.h>

#include "cell/cell.hh"
#include "cell/face.hh"
#include "cell/vertex.hh"

#include "subdiv.hh"
#include "average.hh"

void AvgNOOP::operator()( Cell *cell )
{
	genNormals( cell );
}

void AvgNOOP::genNormal( Vertex *v )
// Generate a normal for Vertex v by averaging together the normals of all of its
// neighboring faces.
{
	Edge *e = v->getEdge();
	Edge *e2 = e->Onext();

	// e,e2 are two successive edges of v.
	Edge *start = e;

	Vec3 vpos = v->pos;
	int cnt = 0;
	Vec3 sumNormal( 0,0,0 );
	do
	{
		++cnt;
					// add into sum the face normal
		sumNormal += norm( cross(e->Dest()->pos-vpos,e2->Dest()->pos-vpos) );
		e = e2;
		e2 = e->Onext();
	} while( e != start );

	v->nor = sumNormal / double(cnt);
}

void AvgNOOP::genNormals( Cell *cell )
// Generate normals for all of the vertices in the cell.
{
	CellVertexIterator verts(cell);
	Vertex *v;
	while( (v = verts.next()) != 0 )
		genNormal( v );
}

void AvgNOOP::applyEvaluation(Cell *cell) {}

void AvgAdHoc::operator()( Cell *cell )
{
	// 1. Generate new positions.  (Put new pos into normal for the time being.)
	{
		CellVertexIterator verts(cell);
		Vertex *v;
		while( (v = verts.next()) != 0 )
			if( !interpolating || v->tag==VODD )
				average( v );
	}
	// 2. Copy positions out of nor into pos.
	{
		CellVertexIterator verts(cell);
		Vertex *v;
		while( (v = verts.next()) != 0 )
			if( !interpolating || v->tag==VODD )
				v->pos = v->nor;
	}
	// 3. Generate new normals.
	AvgNOOP::genNormals( cell );
}

void AvgAdHoc::average( Vertex *v )
// Generate a vertex position for v by averaging the positions of its neighbors.
// Put result in normal field.  (Copy into pos later.)
{
	Edge *start = v->getEdge();
	Edge *e = start;
	v->nor = Vec3(0,0,0);
	int cnt = 0;
	do
	{
		++cnt;
		v->nor += e->Dest()->pos;
		e = e->Onext();
	} while( e != start );

	v->nor /= double(cnt);
}

void AvgAdHoc::applyEvaluation(Cell *cell) {}

/////////////////////////////////////////////////////////////////////////////////////////

void AvgEval::operator()(Cell *cell)
{
	// 1. Generate new positions.  (Put new pos into normal for the time being.)
	{
		CellVertexIterator verts(cell);
		Vertex *v;
		while ((v = verts.next()) != 0)
			if (!interpolating)
				average(v);
			else if (v->tag == VODD)
				butterfly(v);
	}
	// 2. Copy positions out of nor into pos.
	{
		CellVertexIterator verts(cell);
		Vertex *v;
		while ((v = verts.next()) != 0)
			if (!interpolating || v->tag == VODD)
				v->pos = v->nor;
	}
	// 3. Generate new normals.
	AvgNOOP::genNormals(cell); //??????
}

void AvgEval::average(Vertex *v)
// Generate a vertex position for v by averaging the positions of its neighbors.
// Put result in normal field.  (Copy into pos later.)
{
	Edge *start = v->getEdge();
	Edge *e = start;
	v->nor = Vec3(0, 0, 0);
	int n = 0;

	// Iterating through all the vertices
	do
	{
		++n; // number of neighbors
		v->nor += e->Dest()->pos;
		e = e->Onext();
	} while (e != start);

	// Apply Loop Averaging Mask

	// Calculate beta
	double beta = (5.0 / 4.0) - (sqr(3 + 2 * cos((2 * 3.1415) / n))/32);

	// Calculate alpha
	double alpha = (n * (1 - beta)) / beta;

	v->nor += v->pos * alpha;
	// Update Vertex based on averaging mask formula
	v->nor /= double(n + alpha);
}

void AvgEval::applyEvaluation(Cell *cell)
{
	// 1. Generate new positions.  (Put new pos into tmp for the time being.)
	{
		CellVertexIterator verts(cell);
		Vertex *v;
		while ((v = verts.next()) != 0)
			if (!interpolating)
				evaluate(v);
	}
	// 2. Copy positions out of tmp into pos.
	{
		CellVertexIterator verts(cell);
		Vertex *v;
		while ((v = verts.next()) != 0)
			if (!interpolating)
				v->pos = v->tmp;
	}
}

void AvgEval::evaluate(Vertex *v)
// Push a vertex position to the limit using the evaluation mask.
// Put result in tmp field.  (Copy into pos later.)
{
	Edge *start = v->getEdge();
	Edge *e = start;
	v->tmp = Vec3(0, 0, 0);
	int n = 0;

	// Iterating through all the vertices
	do
	{
		++n; // number of neighbors
		v->tmp += e->Dest()->pos;
		e = e->Onext();
	} while (e != start);

	// Apply Loop Evaluation Mask

	// Calculate beta
	double beta = (5.0 / 4.0) - (sqr(3 + 2 * cos((2 * 3.1415) / n)) / 32);

	// Calculate epsilon
	double epsilon = (3 * n) / beta;

	v->tmp += v->pos * epsilon;
	// Update Vertex based on evaluation mask formula
	v->tmp /= double(n + epsilon);
}

void AvgEval::butterfly(Vertex *v)
{
// Generate a vertex position for v by applying averaging masks for butterfly scheme
//    which only takes into account the even vertices
// Put result in normal field.  (Copy into pos later.)
	
	Edge *start = v->getEdge();
	Edge *e = start;
	v->nor = Vec3(0, 0, 0);
	int n = 0;
	Bool flag = true;

	do { // Check for extraordinary neighbors
		n = valence(e->Dest());
		if (!(n == 6))
			flag = false;
		e = e->Onext();
	} while (e != start);
	
	if (flag) // For Regular Interior Point
	{
		// ITERATION LEVEL: 1
		do {
			if (e->Dest()->tag == VEVEN)
				v->nor += 0.5*(e->Dest()->pos);
			e = e->Onext();
		} while (e != start);

		// ITERATION LEVEL: 2
		Edge *start1 = v->getEdge();
		Edge *e1 = start1;
		do {
			Edge *start2 = (e1->Dest())->getEdge();
			Edge *e2 = start2;
			do {
				if (e2->Dest()->tag == VEVEN)
					v->nor += 0.0625*(e2->Dest()->pos);
				e2 = e2->Onext();
			} while (e2 != start2);
			e1 = e1->Onext();
		} while (e1 != start1);

		// ITERATION LEVEL: 3
		Edge *start2 = v->getEdge();
		Edge *e2 = start2;
		do {
			if (e2->Dest()->tag == VEVEN)
			{
				Edge *start3 = (e2->Dest())->getEdge();
				Edge *e3 = start3;

				do {
					e3 = e3->Onext();
				} while (e3->Dest() == v);

				e3 = e3->Onext();
				e3 = e3->Onext();

				{
					Edge *start4 = (e3->Dest())->getEdge();
					Edge *e4 = start4;
					do {
						if (e4->Dest()->tag == VEVEN)
							v->nor += -0.0625*(e3->Dest()->pos);
						e4 = e4->Onext();
					} while (e4 != start4);

				}

				e3 = e3->Onext();
				e3 = e3->Onext();

				{
					Edge *start4 = (e3->Dest())->getEdge();
					Edge *e4 = start4;
					do {
						if (e4->Dest()->tag == VEVEN)
							v->nor += -0.0625*(e3->Dest()->pos);
						e4 = e4->Onext();
					} while (e4 != start4);

				}
			}
			e2 = e2->Onext();
		} while (e2 != start2);
	}

	else {
		do { // Check for extraordinary neighbors
			e = e->Onext();
			n = valence(e->Dest());
		} while (!(n == 6));

		Vertex *v1 = e->Dest();
		Edge *Temp = e->Dest()->getEdge();
		do
		{
			Temp = Temp->Onext();
		} while (Temp->Dest() == v);

		Edge *start1 = Temp;
		Edge *e1 = start1;

		if (n == 3)
		{

			Edge *start2 = e1->Dest()->getEdge();
			Edge *e2 = start2;

			// For s0
			do
			{
				if ((e2->Dest()->tag == VEVEN) && (e2->Dest() != v1))
					v->nor += (5.0 / 12.0)*(e2->Dest()->pos);
			} while (e2 != start2);

			e1 = e1->Onext();

			//For s1 and s2
			do {

				do {
					if ((e2->Dest()->tag == VEVEN) && (e2->Dest() != v1))
						v->nor += (-1.0 / 12.0)*(e2->Dest()->pos);
				} while (e2 != start2);

				e1 = e1->Onext();

			} while (e1 != start1);
		}

		else if (n == 4)
		{
			Edge *start2 = e1->Dest()->getEdge();
			Edge *e2 = start2;

			// For s0			
			do{
				if ((e2->Dest()->tag == VEVEN) && (e2->Dest() != v1))
					v->nor += (3.0 / 8.0)*(e2->Dest()->pos);
			} while (e2 != start2);
			e1 = e1->Onext();
			e1 = e1->Onext();

			//For s2
			do {
				if ((e2->Dest()->tag == VEVEN) && (e2->Dest() != v1))
					v->nor += (-1.0 / 8.0)*(e2->Dest()->pos);
			} while (e2 != start2);
			
		}
		else {
			double j = 0.0;
			do {
				Edge *start2 = e1->Dest()->getEdge();
				Edge *e2 = start2;

				double Sj;
				Sj = 0.25 + cos(2.0*3.1415*j / (double)n) + cos(4.0*j / double(n)) / 2.0;
				Sj /= n;
				do
				{
					if ((e2->Dest()->tag == VEVEN) && (e2->Dest() != v1))
						v->nor += Sj*(e2->Dest()->pos);
				} while (e2 != start2);
				j++;
				e1 = e1->Onext();
			} while (e1 != start1);
		}
	}
		
	//SIR's Method 
	//	          |
	//            |
	//	          V

	/*

	// ITERATION LEVEL: 1
	Edge *start = v->getEdge();
	Edge *e = start;
	v->nor = Vec3(0, 0, 0);
	int n = 0;
	double t = 0.100;
	v->nor += v->pos;

	// Iterating through all the vertices
	do {
		++n; // number of neighbors
			 // v->nor += e->Dest()->pos;
		if (e->Dest()->tag == VODD)
			v->nor += t*(e->Dest()->pos);
		e = e->Onext();
	} while (e != start);

	// ITERATION LEVEL: 2
	Edge *start1 = v->getEdge();
	Edge *e1 = start1;

	do {

		if (e1->Dest()->tag == VEVEN)
		{
			Edge *start2 = (e1->Dest())->getEdge();
			Edge *e2 = start2;

			do {
				e2 = e2->Onext();
			} while (e2->Dest() == v);

			e2 = e2->Onext();
			e2 = e2->Onext();
			v->nor += -t*(e2->Dest()->pos);

			e2 = e2->Onext();
			e2 = e2->Onext();
			v->nor += -t*(e2->Dest()->pos);

		}

		e1 = e1->Onext();

	} while (e1 != start1);
	*/
}
