#include <iostream>     // for cout
#include "Graph2D.h"
#include "cs1037lib-time.h" // for basic timing/pausing
using namespace std;

#define TERMINAL ( (arc *) 1 )		// COPY OF DEFINITION FROM maxflow.cpp
const int INFTY=100; // value of t-links corresponding to hard-constraints 


// function for printing error messages from the max-flow library
inline void mf_assert(char * msg) {cout << "Error from max-flow library: " << msg << endl; cin.ignore(); exit(1);}

// function for computing edge weights from intensity differences
int fn(const double dI, double sigma) {return (int) (5.0/(1.0+(dI*dI/(sigma*sigma))));}  // function used for setting "n-link costs" ("Sensitivity" sigme is a tuning parameter)
const char* image_names2[] = { "logo" , "uwocampus", "canada", "Im_L" , "ImGr_L" , "Michel_L" , "squareCircle_L"}; // an array of image file names for left images
int cameraSpacing[] = {0, 0, 0, 15, 15, 17, 7}; //spacing between each photo, based of user observation. 
bool weightedAverage = false;
bool weightedColours = false;

Graph2D::Graph2D()
: Graph<int,int,int>(0,0,mf_assert)
, m_nodeFlow() 
, m_seeds()
, m_labeling()
{}

Graph2D::Graph2D(Table2D<RGB> & im, double sigma)
: Graph<int,int,int>(im.getWidth()*im.getHeight(),4*im.getWidth()*im.getHeight(),mf_assert)
, m_nodeFlow(im.getWidth(),im.getHeight(),0) 
, m_seeds(im.getWidth(),im.getHeight(),NONE) 
, m_labeling(im.getWidth(),im.getHeight(),NONE)
{
	int wR, wD, height = im.getHeight(), width =  im.getWidth();
	add_node(width*height);    // adding nodes
	for (int y=0; y<(height-1); y++) for (int x=0; x<(width-1); x++) { // adding edges (n-links)
		node_id n = to_node(x,y);
		wR=fn(dI(im[x][y],im[x+1][y]),sigma);   // function dI(...) is declared in Image2D.h
		wD=fn(dI(im[x][y],im[x  ][y+1]),sigma); // function fn(...) is declared at the top of this file
		add_edge( n, n+1,     wR, wR );
		add_edge( n, n+width, wD, wD );
	}
}

/*
* function that when given 2 number it finds the sum between those 2 numbers
*/
int sumOfNums(int num1, int num2) {
	unsigned int i;
	int sum = 0;
	for(i = num1; i <= num2; i++) {
		sum += i;
	}
	return sum;
}

void Graph2D::reset(Table2D<RGB> & im, double sigma, int im_index, int mode, bool eraze_seeds) 
{
	Table2D<RGB> image2 = loadImage<RGB>(to_Cstr(image_names2[im_index] << ".bmp")); // global function defined in Image2D.h
	int camSpace = cameraSpacing[im_index];

	cout << "resetting segmentation" << endl;
	int i, x, y, wR, wD, height = im.getHeight(), width =  im.getWidth();
	double edgeUp = 100000, edgeDown = 100000, tempEdge,  weightedEdgeDown = 0, weightedEdgeUp = 0;
	Point p;
	m_nodeFlow.reset(width,height,0); // NOTE: width and height of m_nodeFlow table 
	                                  // are used in functions to_node() and to_Point(), 
	                                  // Thus, m_nodeFlow MUST BE RESET FIRST!!!!!!!!
	Graph<int,int,int>::reset();

	add_node(width*height);    // adding nodes
	for (y=0; y<(height-1); y++) for (x=0; x<(width-1); x++) { // adding edges (n-links)
		node_id n = to_node(x,y);
		wR=fn(dI(im[x][y],im[x+1][y]),sigma);   // function dI(...) is declared in Image2D.h
		wD=fn(dI(im[x][y],im[x  ][y+1]),sigma); // function fn(...) is declared at the top of this file
		if(mode == 1) {
			wR=150*fn(dI(im[x][y],im[x+1][y]),sigma);   // function dI(...) is declared in Image2D.h
			wD=150*fn(dI(im[x][y],im[x  ][y+1]),sigma); // function fn(...) is declared at the top of this file

			edgeUp = 100000; 
			edgeDown = 100000;

			weightedEdgeUp = 100; 
			weightedEdgeDown = 100;

			int cheapestEdgeDown = 0; //cheapest location of a pixel for a given label 
			int cheapestEdgeUp = 0;

			for(int i = 0; i < camSpace; i++) {
				if((i + x) < image2.getWidth()) {//doesn't go out of bounds
					if(weightedColours) {
						double r = abs(im[x][y].r - image2[x+i][y].r);
						double g = abs(im[x][y].g - image2[x+i][y].g);
						double b = abs(im[x][y].b - image2[x+i][y].b);
						tempEdge = r + g + b;
					} else {
						tempEdge = abs(dI(im[x][y],image2[x+i][y]));
					}
				}
				if(tempEdge < edgeDown) {
					edgeDown = tempEdge;
					cheapestEdgeDown = i;
				}
			}
			for(int i = camSpace; i < (camSpace * 2); i++) {
				if((i + x) < image2.getWidth()) {//doesn't go out of bounds
					if(weightedColours) {
						double r = abs(im[x][y].r - image2[x+i][y].r);
						double g = abs(im[x][y].g - image2[x+i][y].g);
						double b = abs(im[x][y].b - image2[x+i][y].b);
						tempEdge = r + g + b;
					} else {
						tempEdge = abs(dI(im[x][y],image2[x+i][y]));
					}
				}
				if(tempEdge < edgeUp) {
					edgeUp = tempEdge;
					cheapestEdgeUp = i;
				}
			}
			if(weightedAverage) {
				double edgeDownSum = 0; //for weighted arithmetic division
				double edgeUpSum = 0;
				for(int i = 0; i < camSpace; i++) {
					if((i + x) < image2.getWidth()) {//doesn't go out of bounds
						double weight = abs(i - cheapestEdgeDown);
						weightedEdgeDown += (abs(dI(im[x][y],image2[x+i][y])) * weight); //multiply by the weight, which is the distance to the lowest value pixel in the label
						edgeDownSum += weight;
					}
				}
				for(int i = camSpace; i < (camSpace * 2); i++) {
					if((i + x) < image2.getWidth()) {//doesn't go out of bounds
						double weight = abs(i - cheapestEdgeUp);
						weightedEdgeUp += (abs(dI(im[x][y],image2[x+i][y])) * weight);
						edgeUpSum += weight;
					}
				}
				add_tweights(Point(x, y), weightedEdgeDown / edgeDownSum, weightedEdgeUp / edgeUpSum);
			} else {
				add_tweights(Point(x, y), edgeDown, edgeUp);
			}
		}
		add_edge( n, n+1,     wR, wR );
		add_edge( n, n+width, wD, wD );
	}

	m_labeling.reset(width,height,NONE);
	if (eraze_seeds) m_seeds.reset(width,height,NONE); // remove hard constraints
	else { // resets hard "constraints" (seeds)
		cout << "keeping seeds" << endl;
		for (p.y=0; p.y<height; p.y++) for (p.x=0; p.x<width; p.x++) { 
			if      (m_seeds[p]==OBJ) add_tweights( p, INFTY, 0   );
			else if (m_seeds[p]==BKG) add_tweights( p,   0, INFTY );
		}
	}
}

// GUI calls this function when a user left- or right-clicks on an image pixel while in "SEGMENT" mode
void Graph2D::addSeed(Point& p, Label seed_type) 
{
	if (!m_seeds.pointIn(p)) return;
	Label current_constraint = m_seeds[p];
	if (current_constraint==seed_type) return;

	if (current_constraint==NONE) {
		if (seed_type==OBJ) add_tweights( p, INFTY, 0   );
		else                add_tweights( p,   0, INFTY );
	}
	else if (current_constraint==OBJ) {
		if (seed_type==BKG) add_tweights( p, -INFTY, INFTY );
		else                add_tweights( p, -INFTY,   0   );
	}
	else if (current_constraint==BKG) {
		if (seed_type==OBJ) add_tweights( p, INFTY, -INFTY );
		else                add_tweights( p,   0,   -INFTY );
	}
	m_seeds[p]=seed_type;
}

inline Label Graph2D::what_label(Point& p)
{
	node_id i = to_node(p);
	if (nodes[i].parent) return (nodes[i].is_sink) ? BKG : OBJ;
	else                 return NONE;
}

void Graph2D::setLabeling()
{	
	int width  = m_labeling.getWidth();
	int height = m_labeling.getHeight();
	Point p;
	for (p.y=0; p.y<height; p.y++) for (p.x=0; p.x<width; p.x++) m_labeling[p]=what_label(p);
}

int Graph2D::compute_mincut(void (*draw_function)())
{
	int flow = maxflow(draw_function);
	// if visualization of flow is not needed, one can completely remove functions
	// maxflow() and augment() from class Graph2D (t.e. remove their declarations 
	// in Graph2D.h and their implementation code below in Graph2D.cpp) and replace 
	// one line of code above with the following call to the base-class maxflow() method 
	//	      int flow = Graph<int,int,int>::maxflow();

	setLabeling();
	return flow;
}


// overwrites the base class function to allow visualization
// (e.g. calls specified "draw_function" to show each flow augmentation)
int Graph2D::maxflow( void (*draw_function)() ) 
{
	node *i, *j, *current_node = NULL;
	arc *a;
	nodeptr *np, *np_next;

	if (!nodeptr_block)
	{
		nodeptr_block = new DBlock<nodeptr>(NODEPTR_BLOCK_SIZE, error_function);
	}

	changed_list = NULL;
	maxflow_init();

	// main loop
	while ( 1 )
	{
		// test_consistency(current_node);

		if ((i=current_node))
		{
			i -> next = NULL; /* remove active flag */
			if (!i->parent) i = NULL;
		}
		if (!i)
		{
			if (!(i = next_active())) break;
		}

		/* growth */
		if (!i->is_sink)
		{
			/* grow source tree */
			for (a=i->first; a; a=a->next)
			if (a->r_cap)
			{
				j = a -> head;
				if (!j->parent)
				{
					j -> is_sink = 0;
					j -> parent = a -> sister;
					j -> TS = i -> TS;
					j -> DIST = i -> DIST + 1;
					set_active(j);
					add_to_changed_list(j);
				}
				else if (j->is_sink) break;
				else if (j->TS <= i->TS &&
				         j->DIST > i->DIST)
				{
					/* heuristic - trying to make the distance from j to the source shorter */
					j -> parent = a -> sister;
					j -> TS = i -> TS;
					j -> DIST = i -> DIST + 1;
				}
			}
		}
		else
		{
			/* grow sink tree */
			for (a=i->first; a; a=a->next)
			if (a->sister->r_cap)
			{
				j = a -> head;
				if (!j->parent)
				{
					j -> is_sink = 1;
					j -> parent = a -> sister;
					j -> TS = i -> TS;
					j -> DIST = i -> DIST + 1;
					set_active(j);
					add_to_changed_list(j);
				}
				else if (!j->is_sink) { a = a -> sister; break; }
				else if (j->TS <= i->TS &&
				         j->DIST > i->DIST)
				{
					/* heuristic - trying to make the distance from j to the sink shorter */
					j -> parent = a -> sister;
					j -> TS = i -> TS;
					j -> DIST = i -> DIST + 1;
				}
			}
		}

		TIME ++;

		if (a)
		{
			i -> next = i; /* set active flag */
			current_node = i;

			// the line below can be uncommented to visualize current SOURCE and SINK trees
			//... if (draw_function) {setLabeling(); draw_function(); Pause(1);}

			/* augmentation */
			augment(a); // uses overwritten function of the base class (see code below)
			            // that keeps track of flow in table m_nodeFlow
			if (draw_function) {draw_function(); Pause(1);}
			/* augmentation end */

			/* adoption */
			while ((np=orphan_first))
			{
				np_next = np -> next;
				np -> next = NULL;

				while ((np=orphan_first))
				{
					orphan_first = np -> next;
					i = np -> ptr;
					nodeptr_block -> Delete(np);
					if (!orphan_first) orphan_last = NULL;
					if (i->is_sink) process_sink_orphan(i);
					else            process_source_orphan(i);
				}

				orphan_first = np_next;
			}
			/* adoption end */
		}
		else current_node = NULL;
	}
	// test_consistency();

	delete nodeptr_block; 
	nodeptr_block = NULL; 

	maxflow_iteration ++;
	return flow;
}

// overwrites the base class function (to save the nodeFlow in a 2D array m_nodeFlow)
void Graph2D::augment(arc *middle_arc)
{
	node *i;
	arc *a;
	int bottleneck;


	/* 1. Finding bottleneck capacity */
	/* 1a - the source tree */
	bottleneck = middle_arc -> r_cap;
	for (i=middle_arc->sister->head; ; i=a->head)
	{
		a = i -> parent;
		if (a == TERMINAL) break;
		if (bottleneck > a->sister->r_cap) bottleneck = a -> sister -> r_cap;
	}
	if (bottleneck > i->tr_cap) bottleneck = i -> tr_cap;
	/* 1b - the sink tree */
	for (i=middle_arc->head; ; i=a->head)
	{
		a = i -> parent;
		if (a == TERMINAL) break;
		if (bottleneck > a->r_cap) bottleneck = a -> r_cap;
	}
	if (bottleneck > - i->tr_cap) bottleneck = - i -> tr_cap;


	/* 2. Augmenting */
	/* 2a - the source tree */
	middle_arc -> sister -> r_cap += bottleneck;
	middle_arc -> r_cap -= bottleneck;
	for (i=middle_arc->sister->head; ; i=a->head)
	{
		m_nodeFlow[to_point((node_id)(i-nodes))] = 1;
		a = i -> parent;
		if (a == TERMINAL) break;
		a -> r_cap += bottleneck;
		a -> sister -> r_cap -= bottleneck;
		if (!a->sister->r_cap)
		{
			set_orphan_front(i); // add i to the beginning of the adoption list
		}
	}
	i -> tr_cap -= bottleneck;
	if (!i->tr_cap)
	{
		set_orphan_front(i); // add i to the beginning of the adoption list
	}
	/* 2b - the sink tree */
	for (i=middle_arc->head; ; i=a->head)
	{
		m_nodeFlow[to_point((node_id)(i-nodes))] = 1;
		a = i -> parent;
		if (a == TERMINAL) break;
		a -> sister -> r_cap += bottleneck;
		a -> r_cap -= bottleneck;
		if (!a->r_cap)
		{
			set_orphan_front(i); // add i to the beginning of the adoption list
		}
	}
	i -> tr_cap += bottleneck;
	if (!i->tr_cap)
	{
		set_orphan_front(i); // add i to the beginning of the adoption list
	}


	flow += bottleneck;
}
