package layout.algo;

import com.yworks.yfiles.geometry.PointD;
import com.yworks.yfiles.graph.IEdge;
import com.yworks.yfiles.graph.IGraph;
import com.yworks.yfiles.graph.INode;

/**
 * Created by laipple on 05/11/20.
 */
public class SimulatedAnnealingFactory {

    /**
     * Calculate spring forces with classical spring embedder algorithm.
     * @param graph - the input graph.
     * @param lambdaOne - the importance of the node distribution
     * @param lambdaTwo - the importance of the borderline positions
     * @param lambdaThree - the importance of average edge lengths
     * @param lambdaFour - the importance of distances between nodes and edges
     * @param map - the NodeMap, where the calculated forces will be stored (might be non-empty).
     */
    
    static int temperature;
    static int index = 0;
    static double energyOld = 0;
	static double boundThreshold = 10;
	static double bound_top = Double.POSITIVE_INFINITY;
	static double bound_bottom = Double.NEGATIVE_INFINITY;
	static double bound_left = Double.POSITIVE_INFINITY;
	static double bound_right = Double.NEGATIVE_INFINITY;
	static PointD nodePositions[];
    
    public static void simulatedAnnealing (IGraph graph, double lambdaOne, double lambdaTwo, double lambdaThree, double lambdaFour, int iterations)
    {
    	temperature = iterations - index;
    	double positionRadius;
    	
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
    		
    		double dynamic_bound_bottom = bound_top + (40 * graph.getNodes().size());
    		double dynamic_bound_right = bound_left + (40 * graph.getNodes().size());
    		if (dynamic_bound_bottom > bound_bottom) bound_bottom = dynamic_bound_bottom;
    		if (dynamic_bound_right > bound_right) bound_right = dynamic_bound_right;
    	}

    	energyOld = calculateEnergyFunction(graph, lambdaOne, lambdaTwo, lambdaThree, lambdaFour, nodePositions);


    	
    	for (INode n : graph.getNodes())
    	{   
    		double energyNew = 0;   		
    		PointD n_old = nodePositions[(int)n.getTag()];
        	
    		positionRadius = 100 * ((double)temperature / (double)iterations);
    		
        	//Creating randomized coordinates for the new position within 200 units distance in any direction
        	int signx = (Math.random() > 0.5) ? -1 : 1;
        	int signy = (Math.random() > 0.5) ? -1 : 1;
        	double newposx = Math.random() * positionRadius * signx;
        	double newposy = Math.random() * positionRadius * signy;
        	PointD n_new = new PointD((n.getLayout().getCenter().x + newposx), (n.getLayout().getCenter().y + newposy));
        	
        	if(n_new.x > bound_right - boundThreshold || n_new.x < bound_left + boundThreshold || n_new.y > bound_bottom - boundThreshold || n_new.y < bound_top + boundThreshold) continue;
        	
        	nodePositions[(int)n.getTag()] = n_new;
        	
        	/// Now the same procedure is executed for the new Point     
        	energyNew = calculateEnergyFunction(graph, lambdaOne, lambdaTwo, lambdaThree, lambdaFour, nodePositions);
        	
        	if (energyNew < energyOld)
        	{
        		energyOld = energyNew;
            	graph.setNodeCenter(n, n_new);
        	}
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
    	energy += calculateNodeNodeDistances(graph, lambdaOne);
    	energy += calculateBorderlinePositions(graph, lambdaTwo);
    	energy += calculateAvgEdgeLength(graph, lambdaThree);
    	energy += calculateNodeEdgeDistances(graph, lambdaFour);
    	return energy;
    }
    
    public static double calculateNodeNodeDistances (IGraph graph, double lambdaOne)
    {
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
    	return energyNodeDistances * lambdaOne;   	
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
    	/*IListEnumerable<IPort> nPorts = n.getPorts();
		for (IPort p : nPorts)
		{
			for (IEdge e_adj : graph.edgesAt(p, AdjacencyTypes.ALL))
			{
				PointD e_adj_s = new PointD(e_adj.getSourceNode().getLayout().getCenter().x, e_adj.getSourceNode().getLayout().getCenter().y);
				PointD e_adj_t = new PointD(e_adj.getTargetNode().getLayout().getCenter().x, e_adj.getTargetNode().getLayout().getCenter().y);
				cumulativeEdgeDist += e_adj_s.distanceTo(e_adj_t);
				countedges ++;
			}
		}
		if (countedges != 0) averageEdgeDist = cumulativeEdgeDist / countedges;
		
    	cumulativeEdgeDist = 0;
    	countedges = 0;*/
    	
    	double energyEdgeLengths = 0;
    	for (IEdge e : graph.getEdges())
		{			
			energyEdgeLengths += Math.abs(150 - nodePositions[(int)e.getSourceNode().getTag()].distanceTo
					                           (nodePositions[(int)e.getTargetNode().getTag()]));
		}
		energyEdgeLengths /= graph.getEdges().size();
		return energyEdgeLengths * lambdaThree;		
    }
    
    public static double calculateNodeEdgeDistances (IGraph graph, double lambdaFour)
    {
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
		energyNodeEdgeDist /= graph.getNodes().size();
		return energyNodeEdgeDist * lambdaFour;
    }
}
