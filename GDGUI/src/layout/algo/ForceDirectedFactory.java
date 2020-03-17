package layout.algo;

import com.yworks.yfiles.algorithms.GraphConnectivity;
import com.yworks.yfiles.algorithms.YDimension;
import com.yworks.yfiles.algorithms.YOrientedRectangle;
import com.yworks.yfiles.algorithms.YPoint;
import com.yworks.yfiles.algorithms.YVector;
import com.yworks.yfiles.graph.IEdge;
import com.yworks.yfiles.graph.IGraph;
import com.yworks.yfiles.graph.IMapper;
import com.yworks.yfiles.graph.INode;
import com.yworks.yfiles.layout.YGraphAdapter;

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
    	
    	List<YVector> vectors;
    	
    	for (INode u : graph.getNodes())
    	{
    		vectors = new java.util.ArrayList<YVector>();
    		
    		double u_x = u.getLayout().getCenter().x;
            double u_y = u.getLayout().getCenter().y;

            YPoint p_u = new YPoint(u_x, u_y);
            
            for (IEdge z : graph.getEdges())
            {
            	INode z_u = z.getSourceNode();
            	INode z_v = z.getTargetNode();
            	
            	if (u == z_u || u == z_v)
            	{
            		continue;
            	}
            	
            	double z_u_x = z_u.getLayout().getCenter().x;
                double z_u_y = z_u.getLayout().getCenter().y;
                YPoint p_z_u = new YPoint(z_u_x, z_u_y);
                
                double z_v_x = z_v.getLayout().getCenter().x;
                double z_v_y = z_v.getLayout().getCenter().y;
                YPoint p_z_v = new YPoint(z_v_x, z_v_y);
            	
            	YVector z_vec = new YVector(p_z_u, p_z_v);
            	z_vec = YVector.orthoNormal(z_vec);
            	YVector z_vec2 = new YVector(z_vec.rotate(Math.PI));
            	
            	YOrientedRectangle rec1;
            	YOrientedRectangle rec2;
            	
            	YDimension dim = new YDimension(YPoint.distance(p_z_u, p_z_v), YPoint.distance(p_z_u, p_z_v) * 10);
            	
            	if (z_u_x > z_v_x)
            	{
            		if (z_vec.getY() < 0)
                	{
                		rec1 = new YOrientedRectangle(p_z_u, dim, z_vec);
                		rec2 = new YOrientedRectangle(p_z_v, dim, z_vec2);
                	}
            		else
            		{
            			rec1 = new YOrientedRectangle(p_z_v, dim, z_vec);
                		rec2 = new YOrientedRectangle(p_z_u, dim, z_vec2);
            		}
            	}
            	else
            	{
            		if (z_vec.getY() < 0)
            		{
            			rec1 = new YOrientedRectangle(p_z_v, dim, z_vec);
                		rec2 = new YOrientedRectangle(p_z_u, dim, z_vec2);
            		}
            		else
            		{
            			rec1 = new YOrientedRectangle(p_z_u, dim, z_vec);
                		rec2 = new YOrientedRectangle(p_z_v, dim, z_vec2);
            		}
            	}
            	
            	
            	if (rec1.contains(z_u_x, z_u_y) || rec2.contains(z_u_x, z_u_y))
            	{
            		System.out.println("lol");
            	}
            	
            	/**  
            	*   Calculate the shortest distance between the node u and the edge z
            	*   This is always the point (p_z) where the distance vector between u and z is perpendicular to z 
            	*   Use the scale method between p_z and p_u to calculate the electric repulsion and apply 
            	*   the forces to the two nodes of edge z and to u.
            	*/
            }
    	}
    }
}
