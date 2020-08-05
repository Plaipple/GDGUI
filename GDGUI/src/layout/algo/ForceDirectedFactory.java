package layout.algo;

import com.yworks.yfiles.algorithms.GraphConnectivity;
import com.yworks.yfiles.algorithms.LineSegment;
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
    
    
    public static void calculateElectricForcesCrossingResolution (IGraph graph, double electricalRepulsion, double threshold, IMapper<INode, List<YVector>> map)
    {
    	int edgeNumb = 0;
    	boolean[][] edgeCrossings = new boolean[graph.getEdges().size()][graph.getEdges().size()];
    	for (IEdge e : graph.getEdges())
    	{
    		e.setTag(edgeNumb);
    		edgeCrossings[edgeNumb][edgeNumb] = true;
    		edgeNumb ++;
    	}
    	
    	for (IEdge e : graph.getEdges())
    	{
    		YPoint e_s = new YPoint(e.getSourceNode().getLayout().getCenter().x, e.getSourceNode().getLayout().getCenter().y);
    		YPoint e_t = new YPoint(e.getTargetNode().getLayout().getCenter().x, e.getTargetNode().getLayout().getCenter().y);
    		
    		LineSegment ls_e = new LineSegment(e_s, e_t);
    		
    		for (IEdge q : graph.getEdges())
    		{
    		  /*if (q.getSourceNode() == e.getSourceNode() ||
    				q.getSourceNode() == e.getTargetNode() ||
    				q.getTargetNode() == e.getSourceNode() ||
    				q.getTargetNode() == e.getTargetNode()) continue;*/
    			
    			if (edgeCrossings[(int)e.getTag()][(int)q.getTag()]) continue;
    			
    			YPoint q_s = new YPoint(q.getSourceNode().getLayout().getCenter().x, q.getSourceNode().getLayout().getCenter().y);
    			YPoint q_t = new YPoint(q.getTargetNode().getLayout().getCenter().x, q.getTargetNode().getLayout().getCenter().y);
    			
    			LineSegment ls_q = new LineSegment(q_s, q_t);
    			
    			if (LineSegment.getIntersection(ls_e, ls_q) != null)
    			{
    				YVector v_e = new YVector(e_s, e_t);
    				YVector v_q = new YVector(q_s, q_t);
    				
    				
    				double angle = Math.min(Math.min(YVector.angle(v_e, v_q), YVector.angle(v_q, v_e)), Math.min(YVector.angle(v_e.rotate(Math.PI), v_q), YVector.angle(v_q, v_e.rotate(Math.PI))));
    				angle = 180*angle/Math.PI;
    				if (angle == 90) continue;
    				double factor = threshold * electricalRepulsion / Math.pow(angle,2);
    				if (factor > 3) factor = 3;
    				
    				YVector origin_e_t = new YVector(v_e);
    				YVector origin_q_t = new YVector(v_q);
    				origin_e_t.norm();
    				origin_q_t.norm();
    				origin_e_t.scale(factor);
    				origin_q_t.scale(factor);
    				YVector origin_e_s = new YVector(origin_e_t.rotate(Math.PI));
    				YVector origin_q_s = new YVector(origin_q_t.rotate(Math.PI));
    				
    				double angle_es_qs = Math.min(YVector.angle(origin_e_s, origin_q_s), 2*Math.PI - YVector.angle(origin_e_s, origin_q_s));
    				double angle_es_qt = Math.min(YVector.angle(origin_e_s, origin_q_t), 2*Math.PI - YVector.angle(origin_e_s, origin_q_t));
    				
    				if (angle_es_qs < angle_es_qt)    						
    				{
    					map.getValue(e.getSourceNode()).add(origin_q_s);
    					map.getValue(e.getTargetNode()).add(origin_q_t);
    					map.getValue(q.getSourceNode()).add(origin_e_s);
    					map.getValue(q.getTargetNode()).add(origin_e_t);
    				}
    				else
    				{
    					map.getValue(e.getSourceNode()).add(origin_q_t);
    					map.getValue(e.getTargetNode()).add(origin_q_s);
    					map.getValue(q.getSourceNode()).add(origin_e_t);
    					map.getValue(q.getTargetNode()).add(origin_e_s);
    				}    				   	
    				
    				
    				
    				/*boolean e_done = false;
    				boolean q_done = false;
    				boolean xZero = false;
    				boolean yZero = false;
    				
    				if (temp_e.getX() == 0)
    				{
    					if (q_s.y < q_t.y)
    					{
    						if (temp_e.getY() < 0)
    						{
    							map.getValue(q.getSourceNode()).add(temp_e_reverse);
    							map.getValue(q.getTargetNode()).add(temp_e);
    						}
    						else
    						{
    							map.getValue(q.getSourceNode()).add(temp_e);
    							map.getValue(q.getTargetNode()).add(temp_e_reverse);
    						}
    						q_done = true;
    					}
    					else
    					{
    						if (temp_e.getY() < 0)
    						{
    							map.getValue(q.getSourceNode()).add(temp_e);
    							map.getValue(q.getTargetNode()).add(temp_e_reverse);
    						}
    						else
    						{
    							map.getValue(q.getSourceNode()).add(temp_e_reverse);
    							map.getValue(q.getTargetNode()).add(temp_e);
    						}
    						q_done = true;
    					}
    					xZero = true;
    				}
    				else if (temp_e.getY() == 0)
    				{
    					if (e_s.x < e_t.x)
    					{
    						if (temp_q.getX() < 0)
    						{
    							map.getValue(q.getSourceNode()).add(temp_e);
    							map.getValue(q.getTargetNode()).add(temp_e_reverse);
    						}
    						else
    						{
    							map.getValue(q.getSourceNode()).add(temp_e_reverse);
    							map.getValue(q.getTargetNode()).add(temp_e);
    						}
    						q_done = true;
    						yZero = true;
    					}
    					else
    					{
    						if (temp_q.getX() < 0)
    						{
    							map.getValue(q.getSourceNode()).add(temp_e_reverse);
    							map.getValue(q.getTargetNode()).add(temp_e);
    						}
    						else
    						{
    							map.getValue(q.getSourceNode()).add(temp_e);
    							map.getValue(q.getTargetNode()).add(temp_e_reverse);
    						}
    						q_done = true;
    					}
    					yZero = true;
    				}
    				else if (temp_q.getX() == 0)
    				{
    					if (e_s.y < e_t.y)
    					{
    						if (temp_q.getY() < 0)
    						{
    							map.getValue(e.getSourceNode()).add(temp_q);
    							map.getValue(e.getTargetNode()).add(temp_q_reverse);
    						}
    						else
    						{
    							map.getValue(e.getSourceNode()).add(temp_q_reverse);
    							map.getValue(e.getTargetNode()).add(temp_q);
    						}
    						e_done = true;
    					}
    					else
    					{
    						if (temp_q.getY() < 0)
    						{
    							map.getValue(e.getSourceNode()).add(temp_q_reverse);
    							map.getValue(e.getTargetNode()).add(temp_q);
    						}
    						else
    						{
    							map.getValue(e.getSourceNode()).add(temp_q);
    							map.getValue(e.getTargetNode()).add(temp_q_reverse);
    						}
    						e_done = true;
    					}
    					xZero = true;
    				}
    				else if (temp_q.getY() == 0)
    				{
    					if (e_s.x < e_t.x)
    					{
    						if (temp_q.getX() < 0)
    						{
    							map.getValue(e.getSourceNode()).add(temp_q_reverse);
    							map.getValue(e.getTargetNode()).add(temp_q);
    						}
    						else
    						{
    							map.getValue(e.getSourceNode()).add(temp_q);
    							map.getValue(e.getTargetNode()).add(temp_q_reverse);
    						}
    						e_done = true;
    					}
    					else
    					{
    						if (temp_q.getX() < 0)
    						{
    							map.getValue(e.getSourceNode()).add(temp_q);
    							map.getValue(e.getTargetNode()).add(temp_q_reverse);
    						}
    						else
    						{
    							map.getValue(e.getSourceNode()).add(temp_q_reverse);
    							map.getValue(e.getTargetNode()).add(temp_q);
    						}
    						e_done = true;
    					}
    					yZero = true;
    				}
    				
    				if (!xZero)
    				{	
    					if (!e_done && e_s.x < e_t.x)
    					{
    						if (temp_q.getX() < 0)
    						{
    							map.getValue(e.getSourceNode()).add(temp_q_reverse);
    							map.getValue(e.getTargetNode()).add(temp_q);
    						}
    						else
    						{
    							map.getValue(e.getSourceNode()).add(temp_q);
    							map.getValue(e.getTargetNode()).add(temp_q_reverse);
    						}
    					}
    					else if (!e_done && e_s.x > e_t.x)
    					{
    						if (temp_q.getX() < 0)
    						{
    							map.getValue(e.getSourceNode()).add(temp_q);
    							map.getValue(e.getTargetNode()).add(temp_q_reverse);
    						}
    						else
    						{
    							map.getValue(e.getSourceNode()).add(temp_q_reverse);
    							map.getValue(e.getTargetNode()).add(temp_q);
    						}
    					}

    					if (!q_done && q_s.x < q_t.x)
    					{
    						if (temp_e.getX() < 0)
    						{
    							map.getValue(q.getSourceNode()).add(temp_e_reverse);
    							map.getValue(q.getTargetNode()).add(temp_e);
    						}
    						else
    						{
    							map.getValue(q.getSourceNode()).add(temp_e);
    							map.getValue(q.getTargetNode()).add(temp_e_reverse);
    						}
    					}
    					else if (!q_done && q_s.x > q_t.x)
    					{
    						if (temp_e.getX() < 0)
    						{
    							map.getValue(q.getSourceNode()).add(temp_e);
    							map.getValue(q.getTargetNode()).add(temp_e_reverse);
    						}
    						else
    						{
    							map.getValue(q.getSourceNode()).add(temp_e_reverse);
    							map.getValue(q.getTargetNode()).add(temp_e);
    						}
    					}
    				}
    				else if (!yZero)
    				{
    					if (!e_done && e_s.y < e_t.y)
    					{
    						if (temp_q.getY() < 0)
    						{
    							map.getValue(e.getSourceNode()).add(temp_q);
    							map.getValue(e.getTargetNode()).add(temp_q_reverse);
    						}
    						else
    						{
    							map.getValue(e.getSourceNode()).add(temp_q_reverse);
    							map.getValue(e.getTargetNode()).add(temp_q);
    						}
    					}
    					else if (!e_done && e_s.y > e_t.y)
    					{
    						if (temp_q.getY() < 0)
    						{
    							map.getValue(e.getSourceNode()).add(temp_q_reverse);
    							map.getValue(e.getTargetNode()).add(temp_q);
    						}
    						else
    						{
    							map.getValue(e.getSourceNode()).add(temp_q);
    							map.getValue(e.getTargetNode()).add(temp_q_reverse);
    						}
    					}

    					if (!q_done && q_s.y < q_t.y)
    					{
    						if (temp_e.getY() < 0)
    						{
    							map.getValue(q.getSourceNode()).add(temp_e_reverse);
    							map.getValue(q.getTargetNode()).add(temp_e);
    						}
    						else
    						{
    							map.getValue(q.getSourceNode()).add(temp_e);
    							map.getValue(q.getTargetNode()).add(temp_e_reverse);
    						}
    					}
    					else if (!q_done && q_s.y > q_t.y)
    					{
    						if (temp_e.getY() < 0)
    						{
    							map.getValue(q.getSourceNode()).add(temp_e);
    							map.getValue(q.getTargetNode()).add(temp_e_reverse);
    						}
    						else
    						{
    							map.getValue(q.getSourceNode()).add(temp_e_reverse);
    							map.getValue(q.getTargetNode()).add(temp_e);
    						}
    					}
    				}*/
    			}
    			edgeCrossings[(int)e.getTag()][(int)q.getTag()] = true;
    			edgeCrossings[(int)q.getTag()][(int)e.getTag()] = true;
    		}
    	}
    }
    
    
    public static void calculateElectricForcesAngularResolution (IGraph graph, double electricalRepulsion, double threshold, IMapper<INode, List<YVector>> map)
    {
    	for (INode n : graph.getNodes())
    	{
            List<IEdge> edgeList = new ArrayList<IEdge>();
            IListEnumerable<IPort> nPorts = n.getPorts();
            for (IPort p : nPorts)
            {
            	for (IEdge e : graph.edgesAt(p, AdjacencyTypes.ALL))
            	{
            		edgeList.add(e);
            	}
            }
            
            if (edgeList.size()<=1) continue;
            edgeList.sort(new util.CyclicEdgeComparator(graph));
            
            for (int i = 0; i < edgeList.size(); i++)
            {
            	IEdge e1 = (IEdge) edgeList.get(i);
            	IEdge e2 = (IEdge) edgeList.get((i+1)%edgeList.size());
            	YPoint e1_origin;
            	YPoint e1_target;
            	YPoint e2_origin;
            	YPoint e2_target;
            	
            	if (e1.getSourceNode() == n)
            	{
            		e1_origin = new YPoint(e1.getSourceNode().getLayout().getCenter().x, e1.getSourceNode().getLayout().getCenter().y);
            		e1_target = new YPoint(e1.getTargetNode().getLayout().getCenter().x, e1.getTargetNode().getLayout().getCenter().y);
            	}
            	else
            	{
            		e1_origin = new YPoint(e1.getTargetNode().getLayout().getCenter().x, e1.getTargetNode().getLayout().getCenter().y);
            		e1_target = new YPoint(e1.getSourceNode().getLayout().getCenter().x, e1.getSourceNode().getLayout().getCenter().y);
            	}
            	
            	if (e2.getSourceNode() == n)
            	{
                	e2_origin = new YPoint(e2.getSourceNode().getLayout().getCenter().x, e2.getSourceNode().getLayout().getCenter().y);
                	e2_target = new YPoint(e2.getTargetNode().getLayout().getCenter().x, e2.getTargetNode().getLayout().getCenter().y);
            	}
            	else
            	{
            		e2_origin = new YPoint(e2.getTargetNode().getLayout().getCenter().x, e2.getTargetNode().getLayout().getCenter().y);
                	e2_target = new YPoint(e2.getSourceNode().getLayout().getCenter().x, e2.getSourceNode().getLayout().getCenter().y);
            	}
            	
            	YVector vec_e1 = new YVector(e1_target, e1_origin);
            	YVector vec_e2 = new YVector(e2_target, e2_origin);
            	
            	double angle = Math.min(YVector.angle(vec_e1,  vec_e2), 2*Math.PI - YVector.angle(vec_e1, vec_e2));            	
            	angle = 180*angle/Math.PI;
				if (angle == 180) continue;
				double factor = threshold * electricalRepulsion / Math.pow(angle,2);
				if (factor > 3) factor = 3;
				factor *= 2;
				
				vec_e1.norm();
				vec_e2.norm();
				vec_e1.scale(factor);
				vec_e2.scale(factor);
				
				vec_e1 = vec_e1.rotate(Math.PI / 2);
				YVector vec_e1_reverse = new YVector(vec_e1.rotate(Math.PI));
				vec_e2 = vec_e2.rotate(Math.PI / 2);
				YVector vec_e2_reverse = new YVector(vec_e2.rotate(Math.PI));
				
				YPoint potentiale1Node1 = new YPoint((e1_target.x + vec_e1.getX()), (e1_target.y + vec_e1.getY()));
				YPoint potentiale1Node2 = new YPoint((e1_target.x + vec_e1_reverse.getX()), (e1_target.y + vec_e1_reverse.getY()));
				
				YPoint potentiale2Node1 = new YPoint((e2_target.x + vec_e2.getX()), (e2_target.y + vec_e2.getY()));
				YPoint potentiale2Node2 = new YPoint((e2_target.x + vec_e2_reverse.getX()), (e2_target.y + vec_e2_reverse.getY()));
				
				if (potentiale1Node1.distanceTo(e2_target) > potentiale1Node2.distanceTo(e2_target))
				{
					if (e1.getSourceNode() == n)
					{
						map.getValue(e1.getTargetNode()).add(vec_e1);
					}
					else
					{
						map.getValue(e1.getSourceNode()).add(vec_e1);
					}
				}
				else
				{
					if (e1.getSourceNode() == n)
					{
						map.getValue(e1.getTargetNode()).add(vec_e1_reverse);
					}
					else
					{
						map.getValue(e1.getSourceNode()).add(vec_e1_reverse);
					}
				}
				
				if (potentiale2Node1.distanceTo(e1_target) > potentiale2Node2.distanceTo(e1_target))
				{
					if (e2.getSourceNode() == n)
					{
						map.getValue(e2.getTargetNode()).add(vec_e2);
					}
					else
					{
						map.getValue(e2.getSourceNode()).add(vec_e2);
					}
				}
				else
				{
					if (e2.getSourceNode() == n)
					{
						map.getValue(e2.getTargetNode()).add(vec_e2_reverse);
					}
					else
					{
						map.getValue(e2.getSourceNode()).add(vec_e2_reverse);
					}
				}
            }	
    	}
    }
    
    
    public static void angularResolutionForces(IGraph view, double threshold, IMapper<INode, List<YVector>> map )
    {
        double angle;
        List<YVector> vectors;
        for (INode n : view.getNodes()) {
            vectors = new java.util.ArrayList<YVector>();
            //Sort in cyclic order the adjacent edges
            List<IEdge> edgeList = new ArrayList<IEdge>();
            //for (IEdgeCursor e = n.getEdges(); e.ok(); e.next())
            IListEnumerable<IPort> nPorts = n.getPorts();
            for (IPort p : nPorts) {
                for (IEdge e : view.edgesAt(p, AdjacencyTypes.ALL)) {
                    edgeList.add(e);
                }
            }

            if (edgeList.size() == 1) continue;

            edgeList.sort(new util.CyclicEdgeComparator(view));

            for (int i = 0; i < edgeList.size(); i++) {
                //e1=(v,u1) and e2=(v,u2) are adjacent edges at v.

                IEdge e1 = (IEdge) edgeList.get(i);
                IEdge e2 = (IEdge) edgeList.get((i + 1) % edgeList.size());


                INode u1 = (n == e1.getSourceNode()) ? e1.getTargetNode() : e1.getSourceNode();
                INode u2 = (n == e2.getSourceNode()) ? e2.getTargetNode() : e2.getSourceNode();


                YPoint p_v = new YPoint(n.getLayout().getCenter().x, n.getLayout().getCenter().y);
                YPoint p_u1 = new YPoint(u1.getLayout().getCenter().x, u1.getLayout().getCenter().y);
                YPoint p_u2 = new YPoint(u2.getLayout().getCenter().x, u2.getLayout().getCenter().y);
                YPoint p_b = new YPoint(n.getLayout().getCenter().x, n.getLayout().getCenter().y);

                YVector v_u1 = new YVector(p_u1, p_v);
                YVector v_u2 = new YVector(p_u2, p_v);

                // angle is calculated for the rotation of bisector vector force
                angle =  (YVector.angle(v_u2, v_u1)) / 2;
                YVector v_b= v_u2.rotate(angle);

                //Angular force calculation
                YVector u1Force = new YVector(p_u1).orthoNormal(v_b);
                u1Force.norm();
                YVector u2Force = new YVector(p_u2).orthoNormal(v_b);
                u2Force.norm();


                u1Force.scale(threshold *  (-1 / Math.tan(YVector.angle(v_u1, v_b) / 2)));
                map.getValue(u1).add(u1Force);

                u2Force.scale(threshold * (1 / Math.tan(YVector.angle(v_b, v_u2) / 2)));
                map.getValue(u2).add(u2Force);

            }
        }
    }
}
