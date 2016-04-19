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
	AvgNOOP::genNormals(cell);
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
	// if odd selected?
	
	if (v->selected && (evenBetween(v) >= 3))
	{
		v->nor += v->pos;
	}
	
	else if (v->selected && (v->tag == VODD))
	{
		//v->nor += 0.5*v->pos;
		Edge *start1 = v->getEdge();
		Edge *e1 = start1;
		do {
			if (e1->Dest()->selected) {
				
				v->nor += 0.5*e1->Dest()->pos;
			}
			e1 = e1->Onext();
		} while (e1 != start1);
	}

	else if(v->selected && (evenBetween(v)==2))
	{
		v->nor += 0.75*v->pos;
		Edge *start1 = v->getEdge();
		Edge *e1 = start1;
		do {
			if (e1->Dest()->selected)
				v->nor += 0.125*e1->Dest()->pos;
			e1 = e1->Onext();
		} while (e1 != start1);
	}

	

	else {
		do
		{
			++n; // number of neighbors
			v->nor += e->Dest()->pos;
			e = e->Onext();
		} while (e != start);

		// Apply Loop Averaging Mask

		// Calculate beta
		double beta = (5.0 / 4.0) - (sqr(3 + 2 * cos((2 * 3.1415) / n)) / 32);

		// Calculate alpha
		double alpha = (n * (1 - beta)) / beta;

		v->nor += v->pos * alpha;
		// Update Vertex based on averaging mask formula
		v->nor /= double(n + alpha);
	}
}

int AvgEval::evenBetween(Vertex *v)
{
	int count = 0;
	Edge *startTemp = v->getEdge();
	Edge *eTemp = startTemp;
	
	do
	{
		if(eTemp->Dest()->selected == true)
			++count; // number of neighbors
		eTemp = eTemp->Onext();
	} while (eTemp != startTemp);

	return count;
}

void AvgEval::applyEvaluation(Cell *cell)
{
	// 1. Generate normals from tangent masks
	genNormals(cell);

	// 2. Generate new positions.  (Put new pos into tmp for the time being.)
	{
		CellVertexIterator verts(cell);
		Vertex *v;
		while ((v = verts.next()) != 0)
			if (!interpolating)
				evaluate(v);
	}
	// 3. Copy positions out of tmp into pos.
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
	double beta = (5.0 / 4.0) - (sqr(3 + 2 * cos((2 * M_PI) / n)) / 32);

	// Calculate epsilon
	double epsilon = (3 * n) / beta;

	v->tmp += v->pos * epsilon;
	// Update Vertex based on evaluation mask formula
	v->tmp /= double(n + epsilon);
}

void AvgEval::genNormals(Cell *cell)
// Generate normals for all of the vertices in the cell.
{
	CellVertexIterator verts(cell);
	Vertex *v;
	while ((v = verts.next()) != 0)
		genNormal(v);
}

void AvgEval::genNormal(Vertex *v)
// Generate a normal for Vertex v by applying the loop tangent masks.
{
	Edge *start = v->getEdge();
	Edge *e = start;
	Vec3 tan1(0, 0, 0);
	Vec3 tan2(0, 0, 0);
	double n = valence(v);
	double cnt = 0.0;

	// Iterating through all the neighbor vertices
	do
	{
		++cnt; // count of neighbors
		tan1 += tau(cnt, n) * e->Dest()->pos;

		if (cnt == 1)
			tan2 += tau(n, n) * e->Dest()->pos;
		else
			tan2 += tau(cnt - 1, n) * e->Dest()->pos;

		e = e->Onext();
	} while (e != start);

	normalise(tan1);
	normalise(tan2);

	v->nor = cross(tan1, tan2);
	normalise(v->nor);
}

double AvgEval::tau(double i, double n)
{
	return cos((2 * M_PI * i) / n);
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
	int count = 0;

	do { // Check for the extraordinary neighbors
		n = valence(e->Dest());
		if (!(n == 6))
		{
			flag = false;
			count++;
		}
		e = e->Onext();
	} while (e != start);
	
	e = start;

	if (flag) // For a Regular Interior Point
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
				} while (e3->Dest() != v);

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

	else if (count == 1 ){
		do { // For extraordinary neighbors
			e = e->Onext();
			n = valence(e->Dest());
		} while (n == 6);
		Vertex *v1 = e->Dest();
		v->nor += 0.75*v1->pos;
		extraordinary(v, v1, start, e, n);
		v->nor += v->tmp;
	}

	else {		
		do { // For extraordinary neighbors
			e = e->Onext();
			n = valence(e->Dest());
			if (!(n==6))
			{
				Vertex *v1 = e->Dest();
				v->nor += (0.375)*v1->pos;
				extraordinary(v, v1, start, e, n);
				v->nor += (0.5)*(v->tmp);
			}
		} while (e != start);	
	}
}

void AvgEval::extraordinary(Vertex *v, Vertex *v1, Edge *start, Edge *e, int n)
{
	v->tmp = Vec3(0, 0, 0);

	Edge *Temp = e->Dest()->getEdge();
	do
	{
		Temp = Temp->Onext();
	} while (!(Temp->Dest() == v));

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
				v->tmp += (5.0 / 12.0)*(e2->Dest()->pos);
			e2 = e2->Onext();
		} while (e2 != start2);

		e1 = e1->Onext();

		start2 = e1->Dest()->getEdge();
		e2 = start2;

		//For s1
		do {
			if ((e2->Dest()->tag == VEVEN) && (e2->Dest() != v1))
				v->tmp += (-1.0 / 12.0)*(e2->Dest()->pos);
			e2 = e2->Onext();
			} while (e2 != start2);
		e1 = e1->Onext();

		start2 = e1->Dest()->getEdge();
		e2 = start2;

		//For s2
		do {
			if ((e2->Dest()->tag == VEVEN) && (e2->Dest() != v1))
				v->tmp += (-1.0 / 12.0)*(e2->Dest()->pos);
			e2 = e2->Onext();
		} while (e2 != start2);
		e1 = e1->Onext();
		
	}

	else if (n == 4)
	{
		Edge *start2 = e1->Dest()->getEdge();
		Edge *e2 = start2;

		// For s0			
		do {
			if ((e2->Dest()->tag == VEVEN) && (e2->Dest() != v1))
				v->tmp += (3.0 / 8.0)*(e2->Dest()->pos);
			e2 = e2->Onext();
		} while (e2 != start2);


		e1 = e1->Onext();
		e1 = e1->Onext();
		start2 = e1->Dest()->getEdge();
		e2 = start2;

		//For s2
		do {
			if ((e2->Dest()->tag == VEVEN) && (e2->Dest() != v1))
				v->tmp += (-1.0 / 8.0)*(e2->Dest()->pos);
			e2 = e2->Onext();
		} while (e2 != start2);

	}
	else {
		double j = 0.0;
		do {
			Edge *start2 = e1->Dest()->getEdge();
			Edge *e2 = start2;

			double Sj;
			Sj = 0.25 + cos(2.0*3.1415*j / (double)n) + 0.5*cos(4.0*3.1415*j / double(n));
			Sj /= (double)n;
			do
			{
				if ((e2->Dest()->tag == VEVEN) && (e2->Dest() != v1))
					v->tmp += Sj*(e2->Dest()->pos);
				e2 = e2->Onext();
			} while (e2 != start2);
			j++;
			e1 = e1->Onext();
		} while (e1 != start1);
	}
}