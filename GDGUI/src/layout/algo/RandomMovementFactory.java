package layout.algo;

import java.util.ArrayList;
import java.util.List;

import com.yworks.yfiles.algorithms.LineSegment;
import com.yworks.yfiles.algorithms.YPoint;
import com.yworks.yfiles.algorithms.YVector;
import com.yworks.yfiles.geometry.PointD;
import com.yworks.yfiles.graph.AdjacencyTypes;
import com.yworks.yfiles.graph.IEdge;
import com.yworks.yfiles.graph.IGraph;
import com.yworks.yfiles.graph.INode;
import com.yworks.yfiles.graph.IPort;
import com.yworks.yfiles.utils.IListEnumerable;

/**
 * Created by laipple on 05/31/20.
 */

public class RandomMovementFactory 
{	
	/**
     * Calculate new node positions on rays for the Crossing Resolution Algorithm
     * @param graph - the input graph.
     * @param numberRays - the number of rays on which the node can be relocated
     * @param relocateMin - the minimum distance which the new position is set from the old
     * @param relocateMax - the maximum distance which the new position is set from the old
     * @param iterations - the importance of distances between nodes and edges
     */
	
	static double minAngle = Double.POSITIVE_INFINITY;
	static int unchangeTrials = 0;
	static int altModeTrials = 0;
	static double boundThreshold = 20;
	static double bound_top = Double.POSITIVE_INFINITY;
	static double bound_bottom = Double.NEGATIVE_INFINITY;
	static double bound_left = Double.POSITIVE_INFINITY;
	static double bound_right = Double.NEGATIVE_INFINITY;
	static PointD[] nodePositions;
	static boolean[][] edgeCrossings;
	static double[] nodeAngles;
	static double energy = 0;
	static int index = 0;
	
	public static void randomMovement(IGraph graph, int numberRays, int relocateMin, int relocateMax, boolean allNodes, boolean doubleValues, int iterTillAct, int activeIter, int iterations)
	{
		
		if (index == 0)
		{				
			nodePositions = new PointD[graph.getNodes().size()];
			edgeCrossings = new boolean[graph.getEdges().size()][graph.getEdges().size()];
			nodeAngles = new double[graph.getNodes().size()];
			
			nodePositions = initNodeMap(graph, nodePositions);
			edgeCrossings = initEdgeMap(graph, edgeCrossings);
			energy = calculateInitEnergy(graph, edgeCrossings, nodePositions, nodeAngles);
		}	
		
		INode intersectionNodes[] = calculateMinAngleNodes(graph, nodePositions);
		if (intersectionNodes[0] == null)
		{
			index ++;
			if (index >= iterations) index = 0;
			return;
		}
		
		// Alternative Mode if there has been no changes for the last 'iterTillAct' Iterations
		if (unchangeTrials >= iterTillAct)
		{					
			if (allNodes)
			{
				intersectionNodes = enhanceNodeArray(graph);
			}

			if (doubleValues)
			{
				int parameters[] = {numberRays, relocateMin, relocateMax};
				parameters = doubleParameters(parameters);
				numberRays = parameters[0];
				relocateMin = parameters[1];
				relocateMax = parameters[2];
			}
			
			altModeTrials++;
		}
		//----------------------------------------------------------------------------------------

		
		for (int k = 0; k < intersectionNodes.length; k++)
		{
			System.out.println(intersectionNodes[k].toString());
		}
		
		int chosenNode = (int)(Math.random() * intersectionNodes.length);
		INode n = intersectionNodes[chosenNode];
		
		YVector basicVec = new YVector(0, -1);
		double randomRot = Math.random() * (2 * Math.PI);
		basicVec = basicVec.rotate(randomRot);
		double currentAngle = minAngle;
		double localAngle = 0;
		PointD n_new = new PointD();
		PointD n_new_best = new PointD();
		double energyNew = 0;
		double currentEnergy = energy;
				
		
		for (int i = 0; i < numberRays; i++)			
		{
			basicVec.scale((int)(Math.random() * (relocateMax + 1)) + relocateMin);
			basicVec = basicVec.rotate(2 * Math.PI / numberRays);
			n_new = new PointD(basicVec.getX() + n.getLayout().getCenter().x, basicVec.getY() + n.getLayout().getCenter().y);
			
			if(n_new.x > bound_right - boundThreshold || n_new.x < bound_left + boundThreshold || n_new.y > bound_bottom - boundThreshold || n_new.y < bound_top + boundThreshold) continue;
			
			nodePositions[(int)n.getTag()] = n_new;
			
			//localAngle = calculateLowestCrossingAngle(graph, n, n_new);
			energyNew = calculatePositionEnergy(graph, nodePositions, edgeCrossings, nodeAngles, n, n_new);
			if (energyNew > energy)
			{
				n_new_best = n_new;
				energy = energyNew;
			}
			/*
			if (localAngle > minAngle)
			{
				n_new_best = n_new;
				minAngle = localAngle;
			}*/
			
			basicVec.norm();
		
			/*
			basicVec.scale((int)(Math.random() * (relocateMax + 1)) + relocateMin);
			double rotatevalue = 2 * Math.PI / numberRays;
			basicVec = basicVec.rotate(rotatevalue);
			INode testnode = graph.createNode();
			graph.setNodeCenter(testnode, new PointD(basicVec.getX() + n.getLayout().getCenter().x, basicVec.getY() + n.getLayout().getCenter().y));
			basicVec.norm();*/
		}		
		
		if (energy > currentAngle)
		{
			graph.setNodeCenter(n, n_new_best);
			if (altModeTrials >= activeIter)
			{
				unchangeTrials = 0;
				altModeTrials = 0;
			}
			else if (altModeTrials < iterTillAct)
			{
				unchangeTrials = 0;
			}
		}
		else
		{
			nodePositions[(int)n.getTag()] = new PointD(n.getLayout().getCenter().x, n.getLayout().getCenter().y);
			unchangeTrials ++;
			if (altModeTrials >= activeIter)
			{
				unchangeTrials = 0;
				altModeTrials = 0;
			}
		}
		
		
		index ++;
		if (index >= iterations) index = 0;
	}
	
	
	public static double calculatePositionEnergy(IGraph graph, PointD[] nodePositions, boolean[][] edgeCrossings, double[] nodeAngles, INode n, PointD n_p)
	{
		double energy = 0;
		energy += edgeCrossingsEnergy(graph, edgeCrossings);
		
		energy += crossingResolutionEnergy(graph, nodePositions);
		
		INode[] criticalVertices = consideredVerticesNew(graph, n);
		energy += angularResolutionEnergy(graph, nodePositions, criticalVertices, nodeAngles);
		
		return energy;
	}
	
	
	public static double calculateInitEnergy(IGraph graph, boolean[][] edgeCrossings, PointD[] nodePositions, double[] nodeAngles)
	{
		double energy = 0;
		int nodeNumber = 0;
		
		INode[] criticalNodes = new INode[graph.getNodes().size()];
		for (INode n : graph.getNodes())
		{
			criticalNodes[nodeNumber] = n;
			nodeNumber ++;
		}
		
		energy += edgeCrossingsEnergy(graph, edgeCrossings);
		energy += crossingResolutionEnergy(graph, nodePositions);
		energy += angularResolutionEnergy(graph, nodePositions, criticalNodes, nodeAngles);
		return energy;
	}
	
	
	public static INode[] calculateMinAngleNodes(IGraph graph, PointD[] nodePositions)
	{
    	IEdge edges[] = new IEdge[graph.getEdges().size()];
    	INode intersectionNodes[] = new INode[4];
    	double minAngleLocal = Double.POSITIVE_INFINITY;
    	double currentAngle = 0;
    	
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
    				YVector v_n = new YVector(n1, n2);
    	            YVector v_u = new YVector(u1, u2);
    	                	            
    	            currentAngle = Math.min(Math.min(YVector.angle(v_n, v_u), YVector.angle(v_u, v_n)), Math.min(YVector.angle(v_n.rotate(Math.PI), v_u), YVector.angle(v_u, v_n.rotate(Math.PI))));
    	            //System.out.println(180 * currentAngle/Math.PI + "°");
    	            if (currentAngle < minAngleLocal)
    	            {
    	            	minAngleLocal = currentAngle;
    	            	minAngle = currentAngle;
    	            	intersectionNodes[0] = edges[i].getSourceNode();
    	            	intersectionNodes[1] = edges[i].getTargetNode();
    	            	intersectionNodes[2] = edges[j].getSourceNode();
    	            	intersectionNodes[3] = edges[j].getTargetNode();
    	            }
    			}
			}
    	}
    	
    	return intersectionNodes;
	}
	
	
	public static double calculateLowestCrossingAngle (IGraph graph, INode n, PointD n_p)
	{
		double minAngleNew = Double.POSITIVE_INFINITY;
		double currentMinAngleNew = 0;
		IListEnumerable<IPort> nPorts = n.getPorts();
		for (IPort p : nPorts)
		{
			for (IEdge e : graph.edgesAt(p, AdjacencyTypes.ALL))
			{
				YPoint source = new YPoint(n_p.x, n_p.y);
				YPoint n_yp = new YPoint(n.getLayout().getCenter().x, n.getLayout().getCenter().y);
				YPoint target = new YPoint(e.getTargetNode().getLayout().getCenter().x, e.getTargetNode().getLayout().getCenter().y);
				
				if ((target.x == n_yp.x) && (target.y == n_yp.y))
				{
					target = new YPoint(e.getSourceNode().getLayout().getCenter().x, e.getSourceNode().getLayout().getCenter().y);
				}				
				LineSegment l_e1 = new LineSegment(source, target);
				
				for (IEdge q : graph.getEdges())
				{
					if (q.getSourceNode() == n || q.getTargetNode() == n) continue;
					
					YPoint source2 = new YPoint(q.getSourceNode().getLayout().getCenter().x, q.getSourceNode().getLayout().getCenter().y);
				    YPoint target2 = new YPoint(q.getTargetNode().getLayout().getCenter().x, q.getTargetNode().getLayout().getCenter().y);
					LineSegment l_e2 = new LineSegment(source2, target2);
					
					if (LineSegment.getIntersection(l_e2, l_e1) != null)
					{						
						YVector v_n = new YVector(source, target);
	    	            YVector v_u = new YVector(source2, target2);
	    	            
						minAngleNew = Math.min(Math.min(YVector.angle(v_n, v_u), YVector.angle(v_u, v_n)), Math.min(YVector.angle(v_n.rotate(Math.PI), v_u), YVector.angle(v_u, v_n.rotate(Math.PI))));
					}
				}
				
				if (minAngleNew > currentMinAngleNew)
				{
					currentMinAngleNew = minAngleNew;
				}
			}
		}		

		return minAngleNew;
	}
	
	
	public static int[] doubleParameters (int[] parameters)
	{
		for (int i = 0; i < parameters.length; i++)
		{
			parameters[i] *= 2;
		}
		return parameters;
	}
	
	
	public static INode[] enhanceNodeArray (IGraph graph)
	{
		INode[] enhancedNodeArray = new INode[graph.getNodes().size()];
		int i = 0;
		for (INode x : graph.getNodes())
		{
			enhancedNodeArray[i] = x;
			i++;
		}
		return enhancedNodeArray;
	}
	
	
	public static PointD[] initNodeMap (IGraph graph, PointD[] nodePositions)
	{
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
		return nodePositions;
	}
	
	
	public static boolean[][] initEdgeMap (IGraph graph, boolean[][] edgeCrossings)
	{
		int edgeNumber = 0;
		for (IEdge e : graph.getEdges())
		{
			e.setTag(edgeNumber);
			edgeNumber ++;
		}
		
		for (IEdge e : graph.getEdges())
		{
			
			YPoint e_s = new YPoint(e.getSourceNode().getLayout().getCenter().x, e.getSourceNode().getLayout().getCenter().y);
			YPoint e_t = new YPoint(e.getTargetNode().getLayout().getCenter().x, e.getTargetNode().getLayout().getCenter().y);
			
			for (IEdge i : graph.getEdges())
			{
				if (i == e) continue;
				
				YPoint i_s = new YPoint(i.getSourceNode().getLayout().getCenter().x, i.getSourceNode().getLayout().getCenter().y);
				YPoint i_t = new YPoint(i.getTargetNode().getLayout().getCenter().x, i.getTargetNode().getLayout().getCenter().y);
				
				LineSegment lse = new LineSegment (e_s, e_t);
				LineSegment lsi = new LineSegment (i_s, i_t);
				
				if (LineSegment.getIntersection(lse, lsi) != null)
				{
					edgeCrossings[(int)e.getTag()][(int)i.getTag()] = true;
				}
				
			}
		}
		
		return edgeCrossings;
	}
	
	public static INode[] consideredVerticesNew (IGraph graph, INode n)
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

		INode[] consideredVertices = new INode[edgeList.size()+1];
		int k = 0;
		
		for (IEdge e : edgeList)
		{
			consideredVertices[k] = (e.getSourceNode() == n) ? e.getTargetNode() : e.getSourceNode();
			k++;
		}
		consideredVertices[k] = n;

		return consideredVertices;
	}
	
	
	public static boolean[][] edgeCrossingsNew(IGraph graph, boolean[][] edgeCrossings, INode n, PointD nodePositions[])
	{
        IListEnumerable<IPort> nPorts = n.getPorts();
        for (IPort p : nPorts)
        {
        	for (IEdge e : graph.edgesAt(p, AdjacencyTypes.ALL))
        	{
        		YPoint e_s = new YPoint(nodePositions[(int)e.getSourceNode().getTag()].x, nodePositions[(int)e.getSourceNode().getTag()].y);
        		YPoint e_t = new YPoint(nodePositions[(int)e.getTargetNode().getTag()].x, nodePositions[(int)e.getTargetNode().getTag()].y);
        		
        		for (IEdge u : graph.getEdges())
        		{
        			if (u.getSourceNode() == e.getSourceNode() || u.getSourceNode() == e.getTargetNode() || u.getTargetNode() == e.getSourceNode() || u.getTargetNode() == e.getTargetNode()) continue;
        			
        			YPoint u_s = new YPoint(nodePositions[(int)u.getSourceNode().getTag()].x, nodePositions[(int)u.getSourceNode().getTag()].y);
        			YPoint u_t = new YPoint(nodePositions[(int)u.getTargetNode().getTag()].x, nodePositions[(int)u.getTargetNode().getTag()].y);
        			
        			LineSegment lse = new LineSegment (e_s, e_t);
    				LineSegment lsi = new LineSegment (u_s, u_t);
    				
    				if (LineSegment.getIntersection(lse, lsi) != null)
    				{
    					edgeCrossings[(int)e.getTag()][(int)u.getTag()] = true;
    				}
        		}
        	}
        }
                
		return edgeCrossings;
	}
	
	
	public static double edgeCrossingsEnergy(IGraph graph, boolean[][] edgeCrossings)
	{
		double edgeCrossingsEnergy = 0;
		int crossings = 0;
		for (int i = 0; i < edgeCrossings.length; i++)
		{
			for (int j = 0; j < edgeCrossings.length; j++)
			{
				if (edgeCrossings[i][j] == true) crossings ++;
			}
		}
		edgeCrossingsEnergy = (Math.pow(graph.getEdges().size(), 3) / (69 * Math.pow(graph.getNodes().size(), 2))) / crossings;
		return edgeCrossingsEnergy;
	}
	
	
	public static double crossingResolutionEnergy(IGraph graph, PointD[] nodePositions)
	{
		double crossingResolutionEnergy = 0;
		double saveAngle = minAngle;
		calculateMinAngleNodes(graph, nodePositions);
		crossingResolutionEnergy = (180 * minAngle/Math.PI) / 90;
		minAngle = saveAngle;
		return crossingResolutionEnergy;
	}
	
	
	public static double angularResolutionEnergy(IGraph graph, PointD[] nodePositions, INode[] criticalVertices, double[] nodeAngles)
	{
		double angularResolutionEnergy = 0;
		double angle = Double.POSITIVE_INFINITY;
		double currentAngle;
		int edgesCount = 0;
		int nodeIndex = 0;
		
		for (int n = 0; n < criticalVertices.length; n++)
		{
			//Sort in cyclic order the adjacent edges
            List<IEdge> edgeList = new ArrayList<IEdge>();
            IListEnumerable<IPort> nPorts = criticalVertices[n].getPorts();
            for (IPort p : nPorts)
            {
            	for (IEdge e : graph.edgesAt(p, AdjacencyTypes.ALL))
            	{
            		edgeList.add(e);
            	}
            }
            
            currentAngle = Double.POSITIVE_INFINITY;
            
            if (edgeList.size()<=1) continue;            
            edgeList.sort(new util.CyclicEdgeComparator(graph));
            
            for (int i = 0; i < edgeList.size(); i++)
            {
            	IEdge e1 = (IEdge) edgeList.get(i);
            	IEdge e2 = (IEdge) edgeList.get((i+1)%edgeList.size());
            	
            	YPoint e1_s = new YPoint(nodePositions[(int)e1.getSourceNode().getTag()].x, nodePositions[(int)e1.getSourceNode().getTag()].y);
            	YPoint e1_t = new YPoint(nodePositions[(int)e1.getTargetNode().getTag()].x, nodePositions[(int)e1.getTargetNode().getTag()].y);
            	
            	YPoint e2_s = new YPoint(nodePositions[(int)e2.getSourceNode().getTag()].x, nodePositions[(int)e2.getSourceNode().getTag()].y);
            	YPoint e2_t = new YPoint(nodePositions[(int)e2.getTargetNode().getTag()].x, nodePositions[(int)e2.getTargetNode().getTag()].y);

            	YVector v1 = new YVector(e1_s, e1_t);
            	YVector v2 = new YVector(e2_s, e2_t);
            	
            	
            	if (YVector.angle(v1, v2) < currentAngle)
            	{
            		nodeAngles[(int)criticalVertices[n].getTag()] = YVector.angle(v1, v2);
            		currentAngle = YVector.angle(v1, v2);
            	}
            }           
		}
		
		for (int k = 0; k < nodeAngles.length; k++)
		{
			if (nodeAngles[k] < angle && nodeAngles[k] != 0)
			{
				angle = nodeAngles[k];
				nodeIndex = k;
			}
		}
		
		for (INode n : graph.getNodes())
		{
			if ((int)n.getTag() == nodeIndex)
			{
				IListEnumerable<IPort> nPorts = n.getPorts();
				for (IPort p : nPorts)
				{
					for (int l = 0; l < graph.edgesAt(p, AdjacencyTypes.ALL).size(); l++)
					{
						edgesCount++;
					}
				}
				break;
			}
		}
		
		if (angle == Double.POSITIVE_INFINITY) return 0;
		angle = 180 * angle / Math.PI;
		angularResolutionEnergy = angle / (360 / edgesCount);
		return angularResolutionEnergy;
	}
}
