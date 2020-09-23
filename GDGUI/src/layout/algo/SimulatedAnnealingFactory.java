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
     * Different energy functions for the Davidson Harel Simulated Annealing are calculated
     * @param graph - the input graph.
     * @param lambdaOne - the importance of the node distribution
     * @param lambdaTwo - the importance of the borderline positions
     * @param lambdaThree - the importance of average edge lengths
     * @param lambdaFour - the importance of distances between nodes and edges
     * @param area - the area that the graph is allowed to use per node
     */    
    
    public static double calculateNodeNodeDistances (IGraph graph, double lambdaOne, PointD[] nodePositions)
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
    
    
    public static double calculateCrossingNumber (IGraph graph, double lambdaTwo, PointD[] nodePositions)
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
			YPoint n1 = new YPoint(nodePositions[(int)edges[i].getSourceNode().getTag()].x, nodePositions[(int)edges[i].getSourceNode().getTag()].y);
			YPoint n2 = new YPoint(nodePositions[(int)edges[i].getTargetNode().getTag()].x, nodePositions[(int)edges[i].getTargetNode().getTag()].y);
					
    		for (int j = i+1; j < edges.length; j++)
    		{
    			YPoint u1 = new YPoint(nodePositions[(int)edges[j].getSourceNode().getTag()].x, nodePositions[(int)edges[j].getSourceNode().getTag()].y);
    			YPoint u2 = new YPoint(nodePositions[(int)edges[j].getTargetNode().getTag()].x, nodePositions[(int)edges[j].getTargetNode().getTag()].y);
    			
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
    
    
    /*public static double calculateBorderlinePositions (IGraph graph, double lambdaTwo)
    {
    	double energyBorderlines = 0;
    	
    	for (INode u : graph.getNodes())
    	{
    		PointD p_u = nodePositions[(int)u.getTag()];
    		double right = bound_right + boundThreshold - p_u.x;
    		double left = p_u.x - (bound_left - boundThreshold);
    		double top = p_u.y - (bound_top - boundThreshold);
    		double bottom = bound_bottom + boundThreshold - p_u.y;
    		if (right  != 0) right = 1 / Math.pow(right, 2);
    		if (left   != 0) left = 1 / Math.pow(left, 2);
    		if (top    != 0) top = 1 / Math.pow(top, 2);
    		if (bottom != 0) bottom = 1 / Math.pow(bottom, 2);    		
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
    }*/
    
    public static double calculateAvgEdgeLength (IGraph graph, double lambdaThree, PointD[] nodePositions)
    {
    	//calculate the edge lengths of all edges and get the average value
    	
    	double energyEdgeLengths = 0;
    	for (IEdge e : graph.getEdges())
		{			
			energyEdgeLengths += Math.abs(250 - nodePositions[(int)e.getSourceNode().getTag()].distanceTo
					                           (nodePositions[(int)e.getTargetNode().getTag()]));
		}
		energyEdgeLengths /= graph.getEdges().size();
		
		//multiply it by lambdaThree to add the normalizing factor of the result
		return energyEdgeLengths * lambdaThree;		
    }
    
    public static double calculateNodeEdgeDistances (IGraph graph, double lambdaFour, PointD[] nodePositions)
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
