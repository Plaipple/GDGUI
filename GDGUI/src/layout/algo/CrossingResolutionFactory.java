package layout.algo;

import com.yworks.yfiles.algorithms.LineSegment;
import com.yworks.yfiles.algorithms.YPoint;
import com.yworks.yfiles.algorithms.YVector;
import com.yworks.yfiles.geometry.PointD;
import com.yworks.yfiles.graph.IEdge;
import com.yworks.yfiles.graph.IGraph;
import com.yworks.yfiles.graph.INode;

/**
 * Created by laipple on 05/31/20.
 */

public class CrossingResolutionFactory 
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
	
	public static void crossingResolution(IGraph graph, int numberRays, int relocateMin, int relocateMax, boolean allNodes, boolean doubleValues, int iterTillAct, int activeIter, int iterations)
	{
		PointD[] nodePositions = new PointD[graph.getNodes().size()];
		int nodeNumber = 0;
		
		for (INode i : graph.getNodes())
		{   
			i.setTag(nodeNumber);
			nodePositions[(int)i.getTag()] = new PointD(i.getLayout().getCenter().x, i.getLayout().getCenter().y);
			nodeNumber ++;
		}
				
		INode intersectionNodes[] = calculateMinAngleNodes(graph, nodePositions);
		if (intersectionNodes[0] == null) return;
		// Alternative Modus if there has been no changes for the last 'iterTillAct' Iterations
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
		double localAngle = minAngle;
		PointD n_new = new PointD();
		PointD n_new_max = new PointD();
				
		
		for (int i = 0; i < numberRays; i++)			
		{
			basicVec.scale((int)(Math.random() * (relocateMax + 1)) + relocateMin);
			basicVec = basicVec.rotate(2 * Math.PI / numberRays);
			n_new = new PointD(basicVec.getX() + n.getLayout().getCenter().x, basicVec.getY() + n.getLayout().getCenter().y);
			nodePositions[(int)n.getTag()] = n_new;
			
			calculateMinAngleNodes(graph, nodePositions);
			if (minAngle > localAngle)
			{
				n_new_max = n_new;
				localAngle = minAngle;
			}
			basicVec.norm();
		
			/*
			basicVec.scale((int)(Math.random() * (relocateMax + 1)) + relocateMin);
			double rotatevalue = 2 * Math.PI / numberRays;
			basicVec = basicVec.rotate(rotatevalue);
			INode testnode = graph.createNode();
			graph.setNodeCenter(testnode, new PointD(basicVec.getX() + n.getLayout().getCenter().x, basicVec.getY() + n.getLayout().getCenter().y));
			basicVec.norm();*/
		}		
		
		if (localAngle > currentAngle)
		{
			graph.setNodeCenter(n, n_new_max);
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
    	            System.out.println(180 * currentAngle/Math.PI + "°");
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
}
