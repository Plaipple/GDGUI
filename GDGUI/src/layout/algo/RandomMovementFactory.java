package layout.algo;

import java.util.ArrayList;
import java.util.List;

import com.yworks.yfiles.algorithms.Edge;
import com.yworks.yfiles.algorithms.GraphConnectivity;
import com.yworks.yfiles.algorithms.LineSegment;
import com.yworks.yfiles.algorithms.ShortestPaths;
import com.yworks.yfiles.algorithms.YPoint;
import com.yworks.yfiles.algorithms.YVector;
import com.yworks.yfiles.geometry.PointD;
import com.yworks.yfiles.graph.AdjacencyTypes;
import com.yworks.yfiles.graph.IEdge;
import com.yworks.yfiles.graph.IGraph;
import com.yworks.yfiles.graph.INode;
import com.yworks.yfiles.graph.IPort;
import com.yworks.yfiles.layout.YGraphAdapter;
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
	static double[] edgeLengths;
	static double[] nodeAngles;
	static double[][][] shortestPaths;
	static double energy = 0;
	static int index = 0;
	static INode angularNode;
	static INode nodeEdgeNode;
	static double longestEdge = Double.NEGATIVE_INFINITY;
	static double shortestEdge = Double.POSITIVE_INFINITY;
	static INode[] edgeLengthRatioNodes = new INode[4];
	
	public static void randomMovement(IGraph graph, int numberRays, int relocateMin, int relocateMax, boolean allNodes, boolean doubleValues, int iterTillAct, int activeIter, boolean crossResAct, boolean angResAct, boolean numbCrossAct, boolean nodeEdgeResAct, boolean edgeLengthRatioAct, boolean stressAct, int iterations)
	{		
		if (index == 0)
		{
			nodePositions = new PointD[graph.getNodes().size()];
			edgeCrossings = new boolean[graph.getEdges().size()][graph.getEdges().size()];
			nodeAngles = new double[graph.getNodes().size()];
			
			nodePositions = initNodeMap(graph, nodePositions);
			edgeCrossings = initEdgeMap(graph, edgeCrossings);
			if (stressAct) shortestPaths = initShortestPaths(graph, shortestPaths);
			energy = calculateInitEnergy(graph, edgeCrossings, nodePositions, nodeAngles, crossResAct, angResAct, numbCrossAct, nodeEdgeResAct, edgeLengthRatioAct, stressAct);
		}
		
		INode crossingNodes[] = new INode[4];
		if (crossResAct) crossingNodes = calculateMinAngleNodes(graph, nodePositions);
		INode[] intersectionNodes = buildCriticalVerticesArray(graph, crossingNodes, angularNode, nodeEdgeNode, crossResAct, angResAct, nodeEdgeResAct, edgeLengthRatioAct);
		
		
		if (intersectionNodes.length == 0 || intersectionNodes[0] == null)
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

		
	/*	for (int k = 0; k < intersectionNodes.length; k++)
		{
			System.out.println(intersectionNodes[k].toString());
		}
		System.out.println();*/
		
		
		int chosenNode = (int)(Math.random() * intersectionNodes.length);
		INode n = intersectionNodes[chosenNode];
		double positionradiusMax;
		double positionradiusMin;
		
		YVector basicVec = new YVector(0, -1);
		double randomRot = Math.random() * (2 * Math.PI);
		basicVec = basicVec.rotate(randomRot);
		PointD n_new = new PointD();
		PointD n_new_best = new PointD();
		double energyNew = 0;
		double currentEnergy = energy;
		boolean[][] edgeCrossingsSave = edgeCrossings;
		double[] nodeAnglesSave = nodeAngles;
		double[][][] shortestPathsSave = shortestPaths;
		INode[] edgeLengthRatioNodesSave = edgeLengthRatioNodes;
		INode saveAngularNode = angularNode;
		INode saveNodeEdgeNode = nodeEdgeNode;
		INode bestNodeEdgeNode = nodeEdgeNode;
		double longestEdgeSave = longestEdge;
		double shortestEdgeSave = shortestEdge;
		
		for (int i = 0; i < numberRays; i++)			
		{
			positionradiusMax = relocateMax * (iterations - index) / iterations;
			positionradiusMin = relocateMin * (iterations - index) / iterations;
			basicVec.scale((int)(Math.random() * (positionradiusMax + 1)) + positionradiusMin);
			basicVec = basicVec.rotate(2 * Math.PI / numberRays);
			
			n_new = new PointD(basicVec.getX() + n.getLayout().getCenter().x, basicVec.getY() + n.getLayout().getCenter().y);
			
			if(n_new.x > bound_right - boundThreshold || n_new.x < bound_left + boundThreshold || n_new.y > bound_bottom - boundThreshold || n_new.y < bound_top + boundThreshold) continue;
			
			nodePositions[(int)n.getTag()] = n_new;
			
			energyNew = calculatePositionEnergy(graph, nodePositions, edgeCrossings, nodeAngles, shortestPaths, n, n_new, crossResAct, angResAct, numbCrossAct, nodeEdgeResAct, edgeLengthRatioAct, stressAct);
			if (energyNew > energy)
			{
				n_new_best = n_new;
				energy = energyNew;
				nodeEdgeNode = bestNodeEdgeNode;		
			}
			nodePositions[(int)n.getTag()] = new PointD(n.getLayout().getCenter().x, n.getLayout().getCenter().y);
			edgeCrossings = edgeCrossingsSave;
			nodeAngles = nodeAnglesSave;
			angularNode = saveAngularNode;
			nodeEdgeNode = saveNodeEdgeNode;
			edgeLengthRatioNodes = edgeLengthRatioNodesSave;
			shortestEdge = shortestEdgeSave;
			longestEdge = longestEdgeSave;
			shortestPaths = shortestPathsSave;
		
			basicVec.norm();
		
			/*
			basicVec.scale((int)(Math.random() * (relocateMax + 1)) + relocateMin);
			double rotatevalue = 2 * Math.PI / numberRays;
			basicVec = basicVec.rotate(rotatevalue);
			INode testnode = graph.createNode();
			graph.setNodeCenter(testnode, new PointD(basicVec.getX() + n.getLayout().getCenter().x, basicVec.getY() + n.getLayout().getCenter().y));
			basicVec.norm();*/
		}		
		
		if (energy > currentEnergy)
		{
			graph.setNodeCenter(n, n_new_best);
			
			nodePositions[(int)n.getTag()] = n_new_best;
			
			edgeCrossings = edgeCrossingsNew(graph, edgeCrossings, n, nodePositions);
			
			if (stressAct) shortestPaths = shortestPathsNew(graph, shortestPaths, n, nodePositions);
			
			if (edgeLengthRatioAct) edgeLengthRatioEnergy(graph, n);
			
			if (angResAct)
			{
				INode[] criticalVertices = consideredVerticesNew(graph, n);
				angularResolutionEnergy(graph, nodePositions, criticalVertices, nodeAngles);				
			}

			nodeEdgeNode = bestNodeEdgeNode;						
			
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
	
	
	public static double calculatePositionEnergy(IGraph graph, PointD[] nodePositions, boolean[][] edgeCrossings, double[] nodeAngles, double[][][] shortestPaths, INode n, PointD n_p, boolean crossResAct, boolean angResAct, boolean numbCrossAct, boolean nodeEdgeResAct, boolean edgeLengthRatioAct, boolean stressAct)
	{
		double energy = 0;
		PointD old_n_p = new PointD(n.getLayout().getCenter().x, n.getLayout().getCenter().y);
		INode n_copy = graph.createNode(n);
		graph.setNodeCenter(n_copy, n_p);		
		
		if (numbCrossAct)
		{
			edgeCrossings = edgeCrossingsNew(graph, edgeCrossings, n, nodePositions);
			energy += edgeCrossingsEnergy(graph, edgeCrossings);
		}

		
		if (crossResAct) 
			energy += crossingResolutionEnergy(graph, nodePositions);
		
		if (angResAct)
		{
			INode[] criticalVertices = consideredVerticesNew(graph, n);
			energy += angularResolutionEnergy(graph, nodePositions, criticalVertices, nodeAngles);
		}
		
		graph.remove(n_copy);
		if (nodeEdgeResAct)
			energy += nodeEdgeDistanceEnergy(graph);
		
		if (edgeLengthRatioAct)
			energy += edgeLengthRatioEnergy(graph, n);
		
		if (stressAct)
		{
			shortestPaths = shortestPathsNew(graph, shortestPaths, n, nodePositions);
			energy += stressEnergy(graph, shortestPaths);
		}
	
		return energy;
	}
	
	
	public static double calculateInitEnergy(IGraph graph, boolean[][] edgeCrossings, PointD[] nodePositions, double[] nodeAngles, boolean crossResAct, boolean angResAct, boolean numbCrossAct, boolean nodeEdgeResAct, boolean edgeLengthRatioAct, boolean stressAct)
	{
		double energy = 0;
		int nodeNumber = 0;
		
		if (numbCrossAct) energy += edgeCrossingsEnergy(graph, edgeCrossings);
		if (crossResAct) energy += crossingResolutionEnergy(graph, nodePositions);
		
		if (angResAct)
		{
			INode[] criticalNodes = new INode[graph.getNodes().size()];
			for (INode n : graph.getNodes())
			{
				criticalNodes[nodeNumber] = n;
				nodeNumber ++;				
			}
			
			energy += angularResolutionEnergy(graph, nodePositions, criticalNodes, nodeAngles);
		}
		
		if (nodeEdgeResAct) energy += nodeEdgeDistanceEnergy(graph);
		
		if (edgeLengthRatioAct) energy += edgeLengthRatioEnergy(graph, graph.getNodes().first());
		
		if (stressAct) energy += stressEnergy(graph, shortestPaths);
		return energy;
	}
	
	public static INode[] buildCriticalVerticesArray(IGraph graph, INode[] crossingNodes, INode angularNode, INode nodeEdgeNode, boolean crossResAct, boolean angResAct, boolean nodeEdgeResAct, boolean edgeLengthRatioAct)
	{
		int index = 0;
		int arraylength = 0;
		if (crossResAct && crossingNodes[0] != null) arraylength += 4;
		if (angResAct) arraylength += 1;
		if (nodeEdgeResAct) arraylength += 4;
		if (edgeLengthRatioAct) arraylength += 4;
		if (arraylength == 0) arraylength++;
		
		INode criticalVertices[] = new INode[arraylength];
		
		if (crossResAct && crossingNodes[0] != null)
		{
			for (int i = 0; i < crossingNodes.length; i++)
			{
				criticalVertices[i + index] = crossingNodes[i];
			}
			index += 4;
		}
		
		if (angResAct)
		{
			criticalVertices[index] = angularNode;
			index ++;
		}
		
		if (nodeEdgeResAct)
		{
			criticalVertices[index] = nodeEdgeNode;
			criticalVertices[index+1] = nodeEdgeNode;
			criticalVertices[index+2] = nodeEdgeNode;
			criticalVertices[index+3] = nodeEdgeNode;
			index += 4;
		}
		
		if (edgeLengthRatioAct)
		{
			for (int i = 0; i < edgeLengthRatioNodes.length; i++)
			{
				criticalVertices[i + index] = edgeLengthRatioNodes[i];
			}
			index += 4;
		}
		
		if (criticalVertices[0] == null)
		{	
			int random = (int)(Math.random() * graph.getNodes().size());
			
			for (INode n : graph.getNodes())
			{
				if ((int)n.getTag() == random)
				{
					criticalVertices[0] = n;
					break;
				}
				
			}
		}
		
		return criticalVertices;
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
		
		double dynamic_bound_bottom = bound_top + (40 * graph.getNodes().size());
		double dynamic_bound_right = bound_left + (40 * graph.getNodes().size());
		double graph_center;
		
		if (dynamic_bound_bottom > bound_bottom)
		{
			graph_center = bound_top + (bound_bottom - bound_top) / 2;
			bound_bottom = graph_center + (dynamic_bound_bottom - bound_top) / 2;
			bound_top = graph_center - (dynamic_bound_bottom - bound_top) / 2;
		}
		
		if (dynamic_bound_right > bound_right)
		{
			graph_center = bound_left + (bound_right - bound_left) / 2;
			bound_right = graph_center + (dynamic_bound_right - bound_left) / 2;
			bound_left = graph_center - (dynamic_bound_right - bound_left) / 2;
		}
		
		
		return nodePositions;
	}
	
	
	public static boolean[][] initEdgeMap (IGraph graph, boolean[][] edgeCrossings)
	{
		YGraphAdapter adapter = new YGraphAdapter(graph);
		edgeLengths = new double[graph.getEdges().size()];
		int edgeNumber = 0;
		for (IEdge e : graph.getEdges())
		{
			e.setTag(edgeNumber);
			edgeNumber ++;
			
			edgeLengths[adapter.getCopiedEdge(e).index()] = new YPoint(e.getSourceNode().getLayout().getCenter().x, e.getSourceNode().getLayout().getCenter().y).distanceTo
														   (new YPoint(e.getTargetNode().getLayout().getCenter().x, e.getTargetNode().getLayout().getCenter().y));
		}
		
		for (IEdge e : graph.getEdges())
		{			
			YPoint e_s = new YPoint(e.getSourceNode().getLayout().getCenter().x, e.getSourceNode().getLayout().getCenter().y);
			YPoint e_t = new YPoint(e.getTargetNode().getLayout().getCenter().x, e.getTargetNode().getLayout().getCenter().y);
			
			if (YPoint.distance(e_s, e_t) < shortestEdge)
			{
				shortestEdge = YPoint.distance(e_s, e_t);
				edgeLengthRatioNodes[0] = e.getSourceNode();
				edgeLengthRatioNodes[1] = e.getTargetNode();
			}
			
			if (YPoint.distance(e_s, e_t) > longestEdge) 
			{
				longestEdge = YPoint.distance(e_s, e_t);
				edgeLengthRatioNodes[2] = e.getSourceNode();
				edgeLengthRatioNodes[3] = e.getTargetNode();
			}

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
	
	public static double[][][] initShortestPaths (IGraph graph, double[][][] shortestPaths)
	{
		YGraphAdapter adapter = new YGraphAdapter(graph);
		boolean[] reached = new boolean[graph.getNodes().size()];
		shortestPaths = new double[graph.getNodes().size()][graph.getNodes().size()][2];
		
		for (INode n : graph.getNodes())
		{
			YPoint n_yp = new YPoint(nodePositions[(int)n.getTag()].x, nodePositions[(int)n.getTag()].y);
			
			GraphConnectivity.reachable(adapter.getYGraph(), adapter.getCopiedNode(n), false, reached);
			
			for (INode q : graph.getNodes())
			{				
				if (q.getTag() == n.getTag() || !reached[adapter.getCopiedNode(q).index()]) continue;
				
				YPoint q_yp = new YPoint(nodePositions[(int)q.getTag()].x, nodePositions[(int)q.getTag()].y);
				
				if (shortestPaths[(int)n.getTag()][(int)q.getTag()][1] == 0) shortestPaths[(int)n.getTag()][(int)q.getTag()][1] = YPoint.distance(n_yp, q_yp);
				else 														 shortestPaths[(int)q.getTag()][(int)n.getTag()][1] = shortestPaths[(int)n.getTag()][(int)q.getTag()][1];

				if (shortestPaths[(int)n.getTag()][(int)q.getTag()][0] == 0)
				{
			    	Edge[] shortestPathEdges = new Edge[graph.getNodes().size()];
			    	double shortestPathDistance = ShortestPaths.singleSourceSingleSink(adapter.getYGraph(), adapter.getCopiedNode(n), adapter.getCopiedNode(q), false, edgeLengths, shortestPathEdges);
			    	shortestPaths[(int)n.getTag()][(int)q.getTag()][0] = shortestPathDistance;
				}
				else
				{
					shortestPaths[(int)q.getTag()][(int)n.getTag()][0] = shortestPaths[(int)n.getTag()][(int)q.getTag()][0];
				}
			}
		}
		return shortestPaths;
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
    					edgeCrossings[(int)u.getTag()][(int)e.getTag()] = true;
    				}
    				else
    				{
    					edgeCrossings[(int)e.getTag()][(int)u.getTag()] = false;
    					edgeCrossings[(int)u.getTag()][(int)e.getTag()] = false;
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
		crossings /= 2;
		
		if (graph.getEdges().size() > (4 * graph.getNodes().size()))
		{
			edgeCrossingsEnergy = (Math.pow(graph.getEdges().size(), 3) / (69 * Math.pow(graph.getNodes().size(), 2))) / crossings;
		}
		else
		{
			edgeCrossingsEnergy = 1 - (crossings / (Math.pow(graph.getEdges().size(), 2) / 4));
		}
		return edgeCrossingsEnergy;
	}
	
	
	public static double[][][] shortestPathsNew (IGraph graph, double shortestPaths[][][], INode n, PointD nodePositions[])
	{
		YGraphAdapter adapter = new YGraphAdapter(graph);
		boolean[] reached = new boolean[graph.getNodes().size()];
		GraphConnectivity.reachable(adapter.getYGraph(), adapter.getCopiedNode(n), false, reached);
		
		PointD n_save = new PointD(n.getLayout().getCenter().x, n.getLayout().getCenter().y);
		YPoint n_yp = new YPoint(nodePositions[(int)n.getTag()].x, nodePositions[(int)n.getTag()].y);
		graph.setNodeCenter(n, new PointD(n_yp.x, n_yp.y));
		
		for (INode q : graph.getNodes())
		{
			if (q == n || !reached[adapter.getCopiedNode(q).index()]) continue;
			
			YPoint q_yp = new YPoint(nodePositions[(int)q.getTag()].x, nodePositions[(int)q.getTag()].y);
			
			shortestPaths[(int)n.getTag()][(int)q.getTag()][1] = YPoint.distance(n_yp, q_yp);
					
			Edge[] shortestPathEdges = new Edge[graph.getNodes().size()];
			double shortestPathDistance = ShortestPaths.singleSourceSingleSink(adapter.getYGraph(), adapter.getCopiedNode(n), adapter.getCopiedNode(q), false, edgeLengths, shortestPathEdges);
			shortestPaths[(int)n.getTag()][(int)q.getTag()][0] = shortestPathDistance;
		}
		
		graph.setNodeCenter(n, n_save);
		return shortestPaths;
	}
	
	
	public static double stressEnergy (IGraph graph, double shortestPaths[][][])
	{
		double stressEnergy = 0;
		double stressAccu = 0;
		double divisorAccu = 0;
		double maxDist = new YPoint(bound_left, bound_top).distanceTo(new YPoint(bound_right, bound_bottom));
		
		for (int i = 0; i < graph.getNodes().size(); i++)
		{
			for (int j = i+1; j < graph.getNodes().size(); j++)
			{
				if (shortestPaths[i][j][0] == 0) continue;
				stressAccu += Math.pow((1 / Math.pow(shortestPaths[i][j][0], 2)) * shortestPaths[i][j][1] - shortestPaths[i][j][0], 2);
				divisorAccu += Math.pow((1 / Math.pow(shortestPaths[i][j][0], 2)) * maxDist - shortestPaths[i][j][0], 2);
			}
		}
		stressEnergy = 1 - (stressAccu / divisorAccu);
		return stressEnergy;
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

            	YVector v1;
            	YVector v2;
            	
            	if (e1_s.x == e2_s.x && e1_s.y == e2_s.y)
            	{                	
                	v1 = new YVector(e1_s, e1_t);
                	v2 = new YVector(e2_s, e2_t);
            	}
            	else if (e1_s.x == e2_t.x && e1_s.y == e2_t.y)
            	{               	
                	v1 = new YVector(e1_s, e1_t);
                	v2 = new YVector(e2_t, e2_s);
            	}
            	else if (e1_t.x == e2_s.x && e1_t.y == e2_s.y)
            	{
                	v1 = new YVector(e1_t, e1_s);
                	v2 = new YVector(e2_s, e2_t);            		
            	}
            	else
            	{
                	v1 = new YVector(e1_t, e1_s);
                	v2 = new YVector(e2_t, e2_s);            		
            	}
            	
            	
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
				angularNode = n;
				
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
	
	
	public static double nodeEdgeDistanceEnergy(IGraph graph)
    {
		double nodeEdgeDistanceEnergy = 0;
		double currentDistance = 0;
    	double shortestDistance = Double.POSITIVE_INFINITY;
    	double longestDistance = Double.NEGATIVE_INFINITY;
    	double currentShortestDistance = Double.POSITIVE_INFINITY;
    	double currentLongestDistance = Double.NEGATIVE_INFINITY;
    	for (INode n : graph.getNodes())
    	{
    		PointD node = new PointD(nodePositions[(int)n.getTag()].x, nodePositions[(int)n.getTag()].y);
    		for (IEdge e : graph.getEdges())
    		{
    			if (n == e.getSourceNode() || n == e.getTargetNode()) continue;
    			
    			currentDistance = (node.distanceToSegment(new PointD(nodePositions[(int)e.getSourceNode().getTag()].x, nodePositions[(int)e.getSourceNode().getTag()].y), 
    													  new PointD(nodePositions[(int)e.getTargetNode().getTag()].x, nodePositions[(int)e.getTargetNode().getTag()].y)));
    			
    			if (currentDistance < currentShortestDistance) currentShortestDistance = currentDistance;
    			if (currentDistance > longestDistance) longestDistance = currentDistance;
    		}
    		
    		if (currentShortestDistance < shortestDistance)
    		{
    			nodeEdgeNode = n;
    			shortestDistance = currentShortestDistance;
    			//longestDistance = currentLongestDistance;
    			currentShortestDistance = Double.POSITIVE_INFINITY;
    			currentLongestDistance = Double.NEGATIVE_INFINITY;
    		}
    	}
    	
    	nodeEdgeDistanceEnergy = (shortestDistance / longestDistance);
    	return nodeEdgeDistanceEnergy * 5;
    }
	
	
	public static double edgeLengthRatioEnergy (IGraph graph, INode n)
	{
		double edgeLengthRatioEnergy = 0;
		IListEnumerable<IPort> nPorts = n.getPorts();
		for (IPort p : nPorts)
		{
			for (IEdge e : graph.edgesAt(p, AdjacencyTypes.ALL))
			{
				double currentDistance = YPoint.distance(new YPoint(nodePositions[(int)e.getSourceNode().getTag()].x, nodePositions[(int)e.getSourceNode().getTag()].y), new YPoint(nodePositions[(int)e.getTargetNode().getTag()].x, nodePositions[(int)e.getTargetNode().getTag()].y));
				if (currentDistance > shortestEdge && (n == edgeLengthRatioNodes[0] || n == edgeLengthRatioNodes[1]))
				{
					shortestEdge = currentDistance;
					edgeLengthRatioNodes[0] = e.getSourceNode();
					edgeLengthRatioNodes[1] = e.getTargetNode();
				}
				if (currentDistance < longestEdge && (n == edgeLengthRatioNodes[2] || n == edgeLengthRatioNodes[3]))
				{
					longestEdge = currentDistance;
					edgeLengthRatioNodes[2] = e.getSourceNode();
					edgeLengthRatioNodes[3] = e.getTargetNode();
				}
			}
		}

        edgeLengthRatioEnergy = 1 / (longestEdge / shortestEdge);
		return edgeLengthRatioEnergy;
	}
}
