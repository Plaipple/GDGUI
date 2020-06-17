package layout.algo;

import com.yworks.yfiles.algorithms.LineSegment;
import com.yworks.yfiles.algorithms.YPoint;
import com.yworks.yfiles.geometry.PointD;
import com.yworks.yfiles.graph.IEdge;
import com.yworks.yfiles.graph.IGraph;
import com.yworks.yfiles.graph.INode;

/**
 * Created by laipple on 05/11/20.
 */
public class SimulatedAnnealingFactory 
{

    /**
     * Calculate energy function for the Davidson Harel Simulated Annealing
     * @param graph - the input graph.
     * @param lambdaOne - the importance of the node distribution
     * @param lambdaTwo - the importance of the borderline positions
     * @param lambdaThree - the importance of average edge lengths
     * @param lambdaFour - the importance of distances between nodes and edges
     * @param area - the area that the graph is allowed to use per node
     */
    
	//Declaration of static Variables which are also needed to be saved between the iterations
    static int temperature;
    static int index = 0;
    static double energyOld = 0;
	static double boundThreshold = 10;
	static double bound_top = Double.POSITIVE_INFINITY;
	static double bound_bottom = Double.NEGATIVE_INFINITY;
	static double bound_left = Double.POSITIVE_INFINITY;
	static double bound_right = Double.NEGATIVE_INFINITY;
	static PointD nodePositions[];
    
    public static void simulatedAnnealing (IGraph graph, double lambdaOne, double lambdaTwo, double lambdaThree, double lambdaFour, int iterations, int area)
    {
    	//Adjust the temperature value after each iteration
    	temperature = iterations - index;
    	double positionRadius;
    	
    	//if all four lambdas are 0 then none of the energy functions is applied and thus the algorithm can return
    	if (lambdaOne == 0 && lambdaTwo == 0 && lambdaThree == 0 && lambdaFour == 0) return;
    	
    	//Only in the first iteration the bounds are to be set for the area the graph is allowed to use
    	//and the array nodePositions gets filled with positions of the nodes in the graph
    	if (index == 0)
    	{
			nodePositions = new PointD[graph.getNodes().size()];
			int nodeNumber = 0;
						
    		for (INode i : graph.getNodes())
    		{   
    			i.setTag(nodeNumber);
    			bound_top = Math.min(i.getLayout().getCenter().y, bound_top);
    			bound_bottom = Math.max(i.getLayout().getCenter().y, bound_bottom);
    			bound_left = Math.min(i.getLayout().getCenter().x, bound_left);
    			bound_right = Math.max(i.getLayout().getCenter().x, bound_right);
    			nodePositions[(int)i.getTag()] = new PointD(i.getLayout().getCenter().x, i.getLayout().getCenter().y);
    			nodeNumber ++;
    		}
    		
    		//The bottom and right bound are set dynamically depending on how many nodes 
    		//the graph has. The factor 'area' can be set in the panel for the SA Algorithm
    		double dynamic_bound_bottom = bound_top + (area * graph.getNodes().size());
    		double dynamic_bound_right = bound_left + (area * graph.getNodes().size());
    		double graph_center;
    		
    		//Don't take the dynamic values if they are greater than the position of the most outside node
    		//Because it could happen that nodes of the initial graph are placed outside the dynamic bounds
    		//Then these are ignored while running the algorithm and do not change their positions
    		if (dynamic_bound_bottom > bound_bottom)
    		{
    			//calculate the center of the graph's y-axis. Then add half the value of the 
    			//dynamic space for the bottom and subtract half of it for the top
    			graph_center = bound_top + (bound_bottom - bound_top) / 2;
    			bound_bottom = graph_center + (dynamic_bound_bottom - bound_top) / 2;
    			bound_top = graph_center - (dynamic_bound_bottom - bound_top) / 2;
    		}
    		
    		if (dynamic_bound_right > bound_right)
    		{
    			//same as above except it procedures now for the x-axis
    			graph_center = bound_left + (bound_right - bound_left) / 2;
    			bound_right = graph_center + (dynamic_bound_right - bound_left) / 2;
    			bound_left = graph_center - (dynamic_bound_right - bound_left) / 2;
    		}
    		
    		
    	}

    	//Calculate the energy function of the initial graph layout
    	energyOld = calculateEnergyFunction(graph, lambdaOne, lambdaTwo, lambdaThree, lambdaFour, nodePositions);

    	//main loop of the Simulated Annealing algorithm
    	for (INode n : graph.getNodes())
    	{   
    		//Energy function of the new graph positioning after one node has been moved
    		double energyNew = 0;   		
    		PointD n_old = nodePositions[(int)n.getTag()];
        	
    		//The radius of the new position is determined dynamically and decreases during the runtime.
    		//(at 1000 iterations it decreases by 0.1%)
    		positionRadius = 100 * ((double)temperature / (double)iterations);
    		
        	//Creating randomized coordinates for the new position within the distance of positionRadius in any direction
        	int signx = (Math.random() > 0.5) ? -1 : 1;
        	int signy = (Math.random() > 0.5) ? -1 : 1;
        	double newposx = Math.random() * positionRadius * signx;
        	double newposy = Math.random() * positionRadius * signy;
        	PointD n_new = new PointD((n.getLayout().getCenter().x + newposx), (n.getLayout().getCenter().y + newposy));
        	
        	//check if the new position lies within the bounds + some extra space that can be determined manually
        	if(n_new.x > bound_right - boundThreshold || n_new.x < bound_left + boundThreshold || n_new.y > bound_bottom - boundThreshold || n_new.y < bound_top + boundThreshold) continue;
        	
        	//change the position of the node in the array so that the new energy can be calculated for the graph
        	nodePositions[(int)n.getTag()] = n_new;
        	
        	/// Now the same procedure is executed for the new Point
        	energyNew = calculateEnergyFunction(graph, lambdaOne, lambdaTwo, lambdaThree, lambdaFour, nodePositions);
        	
        	//the better positioning leads to saving the coordinates for the node n and updating the actual graph
        	if (energyNew < energyOld)
        	{
        		energyOld = energyNew;
            	graph.setNodeCenter(n, n_new);
        	}
        	//if the new position has a worse energy value the algorithm still might take it
        	//by a probability which is depending on the current temperature and the 
        	//energy difference of the two positions
            else
        	{
        		double ratio = energyOld / energyNew * 0.1;
        	    double temperatureRatio = (double) temperature / (double) iterations;
        		double probability = ratio * temperatureRatio;
        		//System.out.println(probability);
        		if (Math.random() <= probability)
        		{
        			energyOld = energyNew;
                	graph.setNodeCenter(n, n_new);
        		}
        		else
        		{
        			nodePositions[(int)n.getTag()] = n_old;
        		}
        	}    		
    	}
    	
    	index ++;
    	if (temperature <= 1) index = 0;
    }
    
    public static double calculateEnergyFunction(IGraph graph, double lambdaOne, double lambdaTwo, double lambdaThree, double lambdaFour, PointD nodePositions[])
    {
		double energy = 0;		
    	if (lambdaOne != 0) energy += calculateNodeNodeDistances(graph, lambdaOne);
    	//if (lambdaTwo != 0) energy += calculateBorderlinePositions(graph, lambdaTwo);
    	if (lambdaTwo != 0) energy += calculateCrossingNumber(graph, lambdaTwo);
    	if (lambdaThree != 0) energy += calculateAvgEdgeLength(graph, lambdaThree);
    	if (lambdaFour != 0) energy += calculateNodeEdgeDistances(graph, lambdaFour);
    	return energy;
    }
    
    public static double calculateNodeNodeDistances (IGraph graph, double lambdaOne)
    {
    	//calculate the shortest distance to the next node for each node and get the average shortest
    	//distance of all nodes
    	double shortestNodeNodeDist = Double.POSITIVE_INFINITY;
    	double energyNodeDistances = 0;
    	
    	for (INode u : graph.getNodes())
    	{    		   		    		
    		for (INode i : graph.getNodes())
    		{
        		if (i == u || (int)u.getTag() >= (int)i.getTag()) continue;
    			shortestNodeNodeDist = Math.min(shortestNodeNodeDist,
    											nodePositions[(int)u.getTag()].distanceTo(nodePositions[(int)i.getTag()]));
    		}
    		if (shortestNodeNodeDist == Double.POSITIVE_INFINITY) shortestNodeNodeDist = 0;
    		energyNodeDistances += Math.abs(150 - shortestNodeNodeDist);
    		shortestNodeNodeDist = Double.POSITIVE_INFINITY;
    	}
    	energyNodeDistances /= graph.getNodes().size();
    	
    	//add the average value to the energy function multiplied by lambda one which can be chosen
    	//manually in the SA Settings tab
    	return energyNodeDistances * lambdaOne;   	
    }
    
    
    public static double calculateCrossingNumber (IGraph graph, double lambdaTwo)
    {
    	//calculate the number of edge intersections
    	double energyCrossings = 0;
    	int numberOfCrossings = 0;
    	IEdge edges[] = new IEdge[graph.getEdges().size()];
    	
    	int k = 0;
    	for (IEdge e : graph.getEdges())
    	{
    		edges[k] = e;
    		k++;
    	}
    	
    	for (int i = 0; i < edges.length; i++)
    	{
			YPoint n1 = new YPoint(edges[i].getSourceNode().getLayout().getCenter().x, edges[i].getSourceNode().getLayout().getCenter().y);
			YPoint n2 = new YPoint(edges[i].getTargetNode().getLayout().getCenter().x, edges[i].getTargetNode().getLayout().getCenter().y);
					
    		for (int j = i+1; j < edges.length; j++)
    		{
    			YPoint u1 = new YPoint(edges[j].getSourceNode().getLayout().getCenter().x, edges[j].getSourceNode().getLayout().getCenter().y);
    			YPoint u2 = new YPoint(edges[j].getTargetNode().getLayout().getCenter().x, edges[j].getTargetNode().getLayout().getCenter().y);
    			
    			LineSegment l_e1 = new LineSegment(n1, n2);
    			LineSegment l_e2 = new LineSegment(u1, u2);
    			
    			if (LineSegment.getIntersection(l_e2, l_e1) != null)
    			{
    				numberOfCrossings ++;
    			}
    		}
    	}
    	
    	//multiply the result by 5 to get values which are corresponding to the rest of
    	//the energy function operations and make sense
    	//also the lambdaTwo determines how important this operation is respective to the others
    	energyCrossings = numberOfCrossings * 5 * lambdaTwo;    	
    	return energyCrossings;
    }
    
    
    public static double calculateBorderlinePositions (IGraph graph, double lambdaTwo)
    {
    	double energyBorderlines = 0;
    	
    	for (INode u : graph.getNodes())
    	{
    		PointD p_u = nodePositions[(int)u.getTag()];
    		/*double right = bound_right + boundThreshold - p_u.x;
    		double left = p_u.x - (bound_left - boundThreshold);
    		double top = p_u.y - (bound_top - boundThreshold);
    		double bottom = bound_bottom + boundThreshold - p_u.y;
    		if (right  != 0) right = 1 / Math.pow(right, 2);
    		if (left   != 0) left = 1 / Math.pow(left, 2);
    		if (top    != 0) top = 1 / Math.pow(top, 2);
    		if (bottom != 0) bottom = 1 / Math.pow(bottom, 2);*/    		
    		//double formula = right + left + top + bottom;    	
    		double xDistMid = Math.abs(((bound_right + boundThreshold) / 2) - p_u.x);
    		double yDistMid = Math.abs(((bound_bottom + boundThreshold) / 2) - p_u.y);
    		double formula = Math.abs(xDistMid / ((bound_right + boundThreshold) / 2)) + Math.abs(yDistMid / ((bound_bottom + boundThreshold) / 2));
    		//formula = lambdaTwo * formula;
    		//energy += formula;
    		energyBorderlines += formula * 20;
    	}
    	energyBorderlines /= graph.getNodes().size();
    	return energyBorderlines * lambdaTwo;
    }
    
    public static double calculateAvgEdgeLength (IGraph graph, double lambdaThree)
    {
    	//calculate the edge lengths of all edges and get the average value
    	
    	double energyEdgeLengths = 0;
    	for (IEdge e : graph.getEdges())
		{			
			energyEdgeLengths += Math.abs(150 - nodePositions[(int)e.getSourceNode().getTag()].distanceTo
					                           (nodePositions[(int)e.getTargetNode().getTag()]));
		}
		energyEdgeLengths /= graph.getEdges().size();
		
		//multiply it by lambdaThree to add the normalizing factor of the result
		return energyEdgeLengths * lambdaThree;		
    }
    
    public static double calculateNodeEdgeDistances (IGraph graph, double lambdaFour)
    {
    	//For each node calculate the distance between the node and an edge which lies closest to
    	//that node. But the node must not be target- or sourceNode of the edge
    	double shortestNodeEdgeDist = Double.POSITIVE_INFINITY;
    	double energyNodeEdgeDist = 0;
    	
    	for (INode u : graph.getNodes())
		{
			for (IEdge e : graph.getEdges())
    		{
    			if (e.getSourceNode() == u || e.getTargetNode() == u) continue;
    			shortestNodeEdgeDist = Math.min((nodePositions[(int)u.getTag()]).
    											distanceToSegment(nodePositions[(int)e.getSourceNode().getTag()], 
    													          nodePositions[(int)e.getTargetNode().getTag()]),
    											shortestNodeEdgeDist);
    		}
			if (shortestNodeEdgeDist != Double.POSITIVE_INFINITY)
			{
				energyNodeEdgeDist += Math.abs((200 - shortestNodeEdgeDist) * 3);
			}
																	
			shortestNodeEdgeDist = Double.POSITIVE_INFINITY;
		}
    	
    	//calculate the average distance for each node-edge-distance and 
    	//take the lambdaFour into account for the importance of this operation 
    	//respective to the others
		energyNodeEdgeDist /= graph.getNodes().size();
		return energyNodeEdgeDist * lambdaFour;
    }
}
