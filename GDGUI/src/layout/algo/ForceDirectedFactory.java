package layout.algo;

import com.yworks.yfiles.algorithms.GraphConnectivity;
import com.yworks.yfiles.algorithms.YDimension;
import com.yworks.yfiles.algorithms.YOrientedRectangle;
import com.yworks.yfiles.algorithms.YPoint;
import com.yworks.yfiles.algorithms.YVector;
import com.yworks.yfiles.geometry.PointD;
import com.yworks.yfiles.graph.AdjacencyTypes;
import com.yworks.yfiles.graph.IEdge;
import com.yworks.yfiles.graph.IGraph;
import com.yworks.yfiles.graph.IMapper;
import com.yworks.yfiles.graph.INode;
import com.yworks.yfiles.graph.IPort;
import com.yworks.yfiles.layout.YGraphAdapter;
import com.yworks.yfiles.utils.IListEnumerable;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by bekos on 10/28/16.
 */
public class ForceDirectedFactory {

    /**
     * Calculate spring forces with classical spring embedder algorithm.
     * @param graph - the input graph.
     * @param springStiffness - the stiffness of the spring.
     * @param springNaturalLength - the natural length of the spring.
     * @param threshold - a threshold value.
     * @param map - the NodeMap, where the calculated forces will be stored (might be non-empty).
     */
    public static void calculateSpringForcesEades(IGraph graph, double springStiffness, double springNaturalLength, double threshold, IMapper<INode, List<YVector>> map)
    {
        List<YVector> vectors;

        for (INode u : graph.getNodes())
        {
            vectors = new ArrayList<YVector>();

            double u_x = u.getLayout().getCenter().x;
            double u_y = u.getLayout().getCenter().y;

            YPoint p_u = new YPoint(u_x, u_y);

            //Calculate Spring forces...
            for (INode v : graph.neighbors(INode.class, u))
            {
                double v_x = v.getLayout().getCenter().x;
                double v_y = v.getLayout().getCenter().y;

                YPoint p_v = new YPoint(v_x, v_y);

                YVector temp = new YVector(p_v, p_u);
                temp.norm();

                temp.scale(threshold * springStiffness * Math.log(YPoint.distance(p_u, p_v) / springNaturalLength));

                vectors.add(temp);
            }
            map.getValue(u).addAll(vectors);
        }
    }

    /**
     * Calculate electric forces with classical spring embedder algorithm.
     * @param graph - the input graph.
     * @param electricalRepulsion - the electrical repulsion.
     * @param threshold - a threshold value.
     * @param map - the NodeMap, where the calculated forces will be stored (might be non-empty).
     */
    public static void calculateElectricForcesEades(IGraph graph, double electricalRepulsion, double threshold, IMapper<INode, List<YVector>> map)
    {
        YGraphAdapter adapter = new YGraphAdapter(graph);
        boolean[] reached = new boolean[graph.getNodes().size()];

        List<YVector> vectors;

        for (INode u : graph.getNodes())
        {
            vectors = new java.util.ArrayList<YVector>();

            double u_x = u.getLayout().getCenter().x;
            double u_y = u.getLayout().getCenter().y;

            YPoint p_u = new YPoint(u_x, u_y);

            GraphConnectivity.reachable(adapter.getYGraph(), adapter.getCopiedNode(u), false, reached);

            //Calculate Electrical forces
            for (INode v : graph.getNodes())
            {
                if (u == v || !reached[adapter.getCopiedNode(v).index()]) continue;

                double v_x = v.getLayout().getCenter().x;
                double v_y = v.getLayout().getCenter().y;

                YPoint p_v = new YPoint(v_x, v_y);

                YVector temp = new YVector(p_u, p_v);
                temp.norm();

                temp.scale(threshold * electricalRepulsion / Math.pow(YPoint.distance(p_u, p_v),2));
                vectors.add(temp);
            }
            map.getValue(u).addAll(vectors);
        }
    }
    
    /**
     * Calculate electric forces of node to edge with classical spring embedder algorithm.
     * @param graph - the input graph.
     * @param electricalRepulsion - the electrical repulsion.
     * @param threshold - a threshold value.
     * @param map - the NodeMap, where the calculated forces will be stored (might be non-empty).
     */
    public static void calculateElectricForcesNodeEdge(IGraph graph, double electricalRepulsion, double threshold, IMapper<INode, List<YVector>> map)
    {
    	//Declaration of Array Lists for the electric forces that are to be added up, to apply to the several 
    	//target/source nodes of the edges and the nodes which are too close to the edges
    	List<YVector> vectors_node;
    	List<YVector> vectors_edge;
    	
    	//Loop over each edge
    	for (IEdge z : graph.getEdges())
    	{    		
    		vectors_edge = new java.util.ArrayList<YVector>();
    		
    		INode z_u = z.getSourceNode();
    		INode z_v = z.getTargetNode();

    		double z_u_x = z_u.getLayout().getCenter().x;
    		double z_u_y = z_u.getLayout().getCenter().y;
    		YPoint p_z_u = new YPoint(z_u_x, z_u_y);
    		PointD p_z_u_d = new PointD(z_u_x, z_u_y);

    		double z_v_x = z_v.getLayout().getCenter().x;
    		double z_v_y = z_v.getLayout().getCenter().y;
    		YPoint p_z_v = new YPoint(z_v_x, z_v_y);
    		PointD p_z_v_d = new PointD(z_v_x, z_v_y);

    		//Creating the direction-vectors for the rectangles that are used to check, 
    		//if the nodes are in the specific area at which the algorithm should be executed
    		YVector z_vec1 = new YVector(p_z_u, p_z_v);
    		z_vec1 = YVector.orthoNormal(z_vec1);
    		YVector z_vec2 = new YVector(z_vec1.rotate(Math.PI));
    		z_vec2.norm();
    		
    		//Creating the dimension for the two rectangles
    		//the width is always the length of the edge and the height is ten times the edge.    		
    		YDimension dim = new YDimension(YPoint.distance(p_z_u, p_z_v), YPoint.distance(p_z_u, p_z_v) * 10);
        	
    		//The Orientation of the rectangles is that they are always positioned so that the short 
    		//side (width) is the edge and the height is orthogonal to the edge on both sides
        	YOrientedRectangle rec1 = new YOrientedRectangle(p_z_u, dim, z_vec1);
        	YOrientedRectangle rec2 = new YOrientedRectangle(p_z_v, dim, z_vec2);

    		for (INode u : graph.getNodes())
    		{    			
    			//If the iterated node is one of the edges target/source nodes, then 
    			//don't consider that vertex
    			if (u == z_u || u == z_v )
            	{
            		continue;
            	}
    			
    			double factor;
    			vectors_node = new java.util.ArrayList<YVector>();
    			
    			double u_x = u.getLayout().getCenter().x;
                double u_y = u.getLayout().getCenter().y;
    			
                PointD p_u_d = new PointD(u_x, u_y);
                
                //The scale formula for the electric forces is only applied to nodes which
                //are inside one of the two rectangles that were created above
                if (rec1.contains(u_x, u_y) || rec2.contains(u_x, u_y))
            	{            		
            		YVector temp;
            		
            		//Set up the new vector for the force in one of the two directions
            		//relative to the edge depending on which side of the edge the node is positioned
            		if (rec1.contains(u_x, u_y))
            		{
            			temp = new YVector(z_vec1);
            		}
            		else
            		{
            			temp = new YVector(z_vec2);
            		}
            		
            		temp.norm();

            		//Calculate the electrical forces. One for the node that is too close to an edge
            		//And the 180° turned around force for the source/target nodes of the edge
            		//if the factor gets too big create an upper bound, so that the nodes do not jump away
            		//at the beginning of the execution if the graph is very narrow
            		factor = threshold * electricalRepulsion / Math.pow(p_u_d.distanceToSegment(p_z_u_d, p_z_v_d),2);            		
            		if (factor > 3) factor = 3;
            		
                    temp.scale(factor);
                    YVector temp_edge = new YVector(temp.rotate(Math.PI));
                    vectors_node.add(temp);
                    vectors_edge.add(temp_edge);
            		
            	}
                map.getValue(u).addAll(vectors_node);
    		} 
            map.getValue(z_u).addAll(vectors_edge);
            map.getValue(z_v).addAll(vectors_edge);
    	}
    }
    
    
    static int temperature;
    static int index = 0;
    static double energyOld = 0;
	static double boundThreshold = 10;
	static double bound_top = Double.POSITIVE_INFINITY;
	static double bound_bottom = Double.NEGATIVE_INFINITY;
	static double bound_left = Double.POSITIVE_INFINITY;
	static double bound_right = Double.NEGATIVE_INFINITY;
    
    public static void simulatedAnnealing (IGraph graph, double lambdaOne, double lambdaTwo, double lambdaThree, double lambdaFour, int iterations)
    {
    	temperature = iterations - index;
    	
    	if (index == 0)
    	{
        	for (INode i : graph.getNodes())
        	{
        		bound_top = Math.min(i.getLayout().getCenter().y, bound_top);
    			bound_bottom = Math.max(i.getLayout().getCenter().y, bound_bottom);
    			bound_left = Math.min(i.getLayout().getCenter().x, bound_left);
    			bound_right = Math.max(i.getLayout().getCenter().x, bound_right);
        	}
        	
        	energyOld = calculateEnergyFunction(graph, lambdaOne, lambdaTwo, lambdaThree, lambdaFour);

    	}
    	
    	for (INode n : graph.getNodes())
    	{   
    		double energyNew = 0;   		
        	
        	//Creating randomized coordinates for the new position within 200 units distance in any direction
        	int signx = (Math.random() > 0.5) ? -1 : 1;
        	int signy = (Math.random() > 0.5) ? -1 : 1;
        	double newposx = Math.random() * 200 * signx;
        	double newposy = Math.random() * 200 * signy;
        	PointD n_new = new PointD((n.getLayout().getCenter().x + newposx), (n.getLayout().getCenter().y + newposy));
        	
        	if(n_new.x > bound_right - boundThreshold || n_new.x < bound_left + boundThreshold || n_new.y > bound_bottom - boundThreshold || n_new.y < bound_top + boundThreshold) continue;
        	
        	PointD n_old = new PointD(n.getLayout().getCenter().x, n.getLayout().getCenter().y);
        	graph.setNodeCenter(n, n_new);
        	/// Now the same procedure is executed for the new Point     
        	energyNew = calculateEnergyFunction(graph, lambdaOne, lambdaTwo, lambdaThree, lambdaFour);
        	
        	if (energyNew < energyOld)
        	{
        		energyOld = energyNew;
        	}
            else
        	{
        		//double difference = energyOld - energyNew;
        		double ratio = energyOld / energyNew * 0.1;
        		//double probability = Math.exp((energyOld - energyNew)); // (double) temperature);#
        	    double temperatureRatio = (double) temperature / (double) iterations;
        		double probability = ratio * temperatureRatio;
        		//probability = (probability - 1) * -1;
        		System.out.println(probability);
        		if (Math.random() <= probability)
        		{
        			energyOld = energyNew;
        		}
        		else
        		{
        			graph.setNodeCenter(n, n_old);
        		}
        	}    		
    	}
    	
    	index ++;
    	if (temperature <= 0) index = 0;
    }
    
    public static double calculateEnergyFunction(IGraph graph, double lambdaOne, double lambdaTwo, double lambdaThree, double lambdaFour)
    {
		double energy = 0;
    	double shortestNodeEdgeDist = Double.POSITIVE_INFINITY;
    	double shortestNodeNodeDist = Double.POSITIVE_INFINITY;
    	double energyNodeDistances = 0;
    	double energyBorderlines = 0;
    	double energyEdgeLengths = 0;
    	double energyNodeEdgeDist = 0;
    	
    	//Distances between each pair of nodes
    	for (INode u : graph.getNodes())
    	{    		
    		for (INode i : graph.getNodes())
    		{
        		if (i == u) continue;
    			PointD p_u = new PointD(u.getLayout().getCenter().x, u.getLayout().getCenter().y);
    			PointD p_i = new PointD(i.getLayout().getCenter().x, i.getLayout().getCenter().y);
    			
    			//energy += lambdaOne / Math.pow(p_u.distanceTo(p_i), 2);
    			shortestNodeNodeDist = Math.min(shortestNodeNodeDist, p_u.distanceTo(p_i));    			
    		}
    		energyNodeDistances += Math.abs(150 - shortestNodeNodeDist);
    		shortestNodeNodeDist = Double.POSITIVE_INFINITY;
    	}
    	energyNodeDistances /= graph.getNodes().size();
    	energyNodeDistances *= lambdaOne;
    	energy += energyNodeDistances;
		
    	//Borderline Thresholds
    	for (INode u : graph.getNodes())
    	{
    		PointD p_n = new PointD(u.getLayout().getCenter().x, u.getLayout().getCenter().y);
    		double right = bound_right + boundThreshold - p_n.x;
    		double left = p_n.x - (bound_left - boundThreshold);
    		double top = p_n.y - (bound_top - boundThreshold);
    		double bottom = bound_bottom + boundThreshold - p_n.y;
    		/*if (right  != 0) right = 1 / Math.pow(right, 2);
    		if (left   != 0) left = 1 / Math.pow(left, 2);
    		if (top    != 0) top = 1 / Math.pow(top, 2);
    		if (bottom != 0) bottom = 1 / Math.pow(bottom, 2);*/    		
    		//double formula = right + left + top + bottom;    	
    		double xDistMid = Math.abs(((bound_right + boundThreshold) / 2) - p_n.x);
    		double yDistMid = Math.abs(((bound_bottom + boundThreshold) / 2) - p_n.y);
    		double formula = Math.abs(xDistMid / ((bound_right + boundThreshold) / 2)) + Math.abs(yDistMid / ((bound_bottom + boundThreshold) / 2));
    		//formula = lambdaTwo * formula;
    		//energy += formula;
    		energyBorderlines += formula * 20;
    	}
    	energyBorderlines /= graph.getNodes().size();
    	energyBorderlines *= lambdaTwo;
    	//energy += energyBorderlines;
		
		
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
		
    	//Edge Lengths
		for (IEdge e : graph.getEdges())
		{
			PointD e_s = new PointD(e.getSourceNode().getLayout().getCenter().x, e.getSourceNode().getLayout().getCenter().y);
			PointD e_t = new PointD(e.getTargetNode().getLayout().getCenter().x, e.getTargetNode().getLayout().getCenter().y);
			//energy += lambdaThree/(Math.pow(e_s.distanceTo(e_t), 2));
			energyEdgeLengths += Math.abs(150 - e_s.distanceTo(e_t));
		}
		energyEdgeLengths /= graph.getEdges().size();
		energyEdgeLengths *= lambdaThree;
	    energy += energyEdgeLengths;
		
		//Node to Edge Distances
		for (INode u : graph.getNodes())
		{
			for (IEdge e : graph.getEdges())
    		{
    			if (e.getSourceNode() == u || e.getTargetNode() == u) continue;
    			shortestNodeEdgeDist = Math.min(new PointD(u.getLayout().getCenter().x, u.getLayout().getCenter().y).
    											distanceToSegment(new PointD(e.getSourceNode().getLayout().getCenter().x, e.getSourceNode().getLayout().getCenter().y), 
    													          new PointD(e.getTargetNode().getLayout().getCenter().x, e.getTargetNode().getLayout().getCenter().y)),
    											shortestNodeEdgeDist);
    		}
			if (shortestNodeEdgeDist != Double.POSITIVE_INFINITY)
			{
				//energy += lambdaFour / (Math.pow(shortestNodeEdgeDist, 2));
				if (Math.abs(150 - shortestNodeEdgeDist) > 50) 
				{
					energyNodeEdgeDist += Math.abs((150 - shortestNodeEdgeDist) * 3);
				}
				else
				{
					energyNodeEdgeDist += Math.abs(150 - shortestNodeEdgeDist);
				}
			}
																	
			shortestNodeEdgeDist = Double.POSITIVE_INFINITY;
		}
		energyNodeEdgeDist /= graph.getNodes().size();
		energyNodeEdgeDist *= lambdaFour;
		energy += energyNodeEdgeDist;
    	return energy;
    }
}
