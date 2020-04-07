package util;

import java.util.ArrayList;
import java.util.List;

import com.yworks.yfiles.algorithms.Edge;
import com.yworks.yfiles.algorithms.EdgeList;
import com.yworks.yfiles.algorithms.IEdgeCursor;
import com.yworks.yfiles.algorithms.INodeCursor;
import com.yworks.yfiles.algorithms.LineSegment;
import com.yworks.yfiles.algorithms.Node;
import com.yworks.yfiles.algorithms.YPoint;
import com.yworks.yfiles.algorithms.YVector;
import com.yworks.yfiles.geometry.PointD;
import com.yworks.yfiles.graph.AdjacencyTypes;
import com.yworks.yfiles.graph.IEdge;
import com.yworks.yfiles.graph.INode;
import com.yworks.yfiles.graph.IPort;
import com.yworks.yfiles.utils.IListEnumerable;
import com.yworks.yfiles.view.GraphComponent;

/**
 *
 * @author michael
 */
public class Utilities {

    public static void displayGraphRandomly(GraphComponent graph)
    {
        java.util.Random r = new java.util.Random(System.currentTimeMillis());

        //for (INodeCursor nc = graph.getGraph(); nc.ok(); nc.next()) 
        for (INode n : graph.getGraph().getNodes())
        {
        	graph.getGraph().setNodeCenter(n, new PointD(r.nextInt((int)graph.getViewport().getWidth()), r.nextInt((int)graph.getViewport().getHeight())));
        	
            //graph.getGraph().createNode(new PointD(r.nextInt((int)graph.getViewport().getWidth()), r.nextInt((int)graph.getViewport().getHeight())));
        }
    }

    public static void printDistances(GraphComponent graph)
    {
        //for (IEdgeCursor ec = graph.getEdges(); ec.ok(); ec.next())
        for (IEdge e : graph.getGraph().getEdges())
        {
            System.out.print("From " + e.getSourceNode() + " To " + e.getTargetNode() + " : ");
            System.out.println(YPoint.distance(new YPoint(e.getSourceNode().getLayout().getCenter().x, (e.getSourceNode().getLayout().getCenter().y)), new YPoint(e.getTargetNode().getLayout().getCenter().x, (e.getTargetNode().getLayout().getCenter().y))));
        }
        System.out.println();
    }

    public static YVector bisectionVector(YVector a, YVector b)
    {
        YVector unit_a = YVector.getNormal(a);
        YVector unit_b = YVector.getNormal(b);
        YVector bisection = YVector.add(unit_a, unit_b);
        bisection.norm();

        return bisection;
    }

    public static void swap(Object a, Object b)
    {
        Object temp = a;
        a = b;
        b = temp;
    }

    public static int maxDegree(GraphComponent graph)
    {
        int maxDegree = -1;
        int countedges;
        for (INode n : graph.getGraph().getNodes())
        {
        	countedges = 0;
        	IListEnumerable<IPort> nPorts = n.getPorts();
        	for (IPort p : nPorts)
            {
        		countedges += graph.getGraph().edgesAt(p, AdjacencyTypes.ALL).size();
            }
    		maxDegree = Math.max(maxDegree, countedges);
        }
        return maxDegree;
    }

    public static double calculateAngularResolution(GraphComponent view)
    {
        double angle = Double.MAX_VALUE;

        for (INode n : view.getGraph().getNodes())
        {
            //Sort in cyclic order the adjacent edges
            List<IEdge> edgeList = new ArrayList<IEdge>();
            //for (IEdgeCursor e = n.getEdges(); e.ok(); e.next())
            IListEnumerable<IPort> nPorts = n.getPorts();
            for (IPort p : nPorts)
            {
            	for (IEdge e : view.getGraph().edgesAt(p, AdjacencyTypes.ALL))
            	{
            		edgeList.add(e);
            	}
            }

            if (edgeList.size()==1) continue;

            edgeList.sort(new util.CyclicEdgeComparator(view.getGraph()));

            for (int i = 0; i < edgeList.size(); i++)
            {
               //e1=(v,u1) and e2=(v,u2) are adjacent edges at v.

               IEdge e1 = (IEdge) edgeList.get(i);
               IEdge e2 = (IEdge) edgeList.get((i+1)%edgeList.size());


               INode u1 = ( n == e1.getSourceNode() ) ? e1.getTargetNode() : e1.getSourceNode();
               INode u2 = ( n == e2.getSourceNode() ) ? e2.getTargetNode() : e2.getSourceNode();


               YPoint p_v = new YPoint(n.getLayout().getCenter().x, n.getLayout().getCenter().y);
               YPoint p_u1 = new YPoint(u1.getLayout().getCenter().x, u1.getLayout().getCenter().y);
               YPoint p_u2 = new YPoint(u2.getLayout().getCenter().x, u2.getLayout().getCenter().y);

               YVector v_u1 = new YVector(p_u1, p_v);
               YVector v_u2 = new YVector(p_u2, p_v);

               angle = Math.min(angle, YVector.angle(v_u2, v_u1));
            }
        }
        return 180*angle/Math.PI;
    }
    
    public static double calculateAverageEdgeLength(GraphComponent view)
    {
        double totalEdgeLength = 0.0;
        for (IEdge e : view.getGraph().getEdges())
        {   
        	/*
            Node u1 = (e.getSourceNode().getLayout().getCenter().x <= view.getRealizer(e.target()).getCenterX()) ? e.source() : e.target();
            Node u2 = (view.getRealizer(e.source()).getCenterX() >  view.getRealizer(e.target()).getCenterX()) ? e.source() : e.target();
           
            YPoint p_u1 = new YPoint(view.getRealizer(u1).getCenterX(), view.getRealizer(u1).getCenterY());
            YPoint p_u2 = new YPoint(view.getRealizer(u2).getCenterX(), view.getRealizer(u2).getCenterY());
            */

        	INode source = e.getSourceNode();
        	INode target = e.getTargetNode();
        	double centerx = source.getLayout().getCenter().x;
        	double targety = target.getLayout().getCenter().y;
        	
        	LineSegment ls = new LineSegment(new YPoint(e.getSourceNode().getLayout().getCenter().x, e.getSourceNode().getLayout().getCenter().y), 
        								 	 new YPoint(e.getTargetNode().getLayout().getCenter().x, e.getTargetNode().getLayout().getCenter().y));
        	
        	
        	
            //totalEdgeLength += YPoint.distance(p_u1,  p_u2);
        	totalEdgeLength += ls.length();
        }
        return totalEdgeLength / view.getGraph().getEdges().size();
    }

    public static int calculateNumberOfCrossings(GraphComponent view)
    {
        int numberOfCrossings = 0;

        IEdge[] edgeArray = new IEdge[view.getGraph().getEdges().size()];

        int k=0;
        for (IEdge e : view.getGraph().getEdges())
        {
            edgeArray[k] = e;
            k++;
        }
        for (int i=0; i<edgeArray.length; i++)
        {
            for (int j=i+1; j<edgeArray.length; j++)
            {
                // e1 = (u1,u2) u1<u2 and e2 = (v1,v2) v1<v2

                INode u1 = (edgeArray[i].getSourceNode().getLayout().getCenter().x <= edgeArray[i].getTargetNode().getLayout().getCenter().x) ? edgeArray[i].getSourceNode() : edgeArray[i].getTargetNode();
                INode u2 = (edgeArray[i].getSourceNode().getLayout().getCenter().x >  edgeArray[i].getTargetNode().getLayout().getCenter().x) ? edgeArray[i].getSourceNode() : edgeArray[i].getTargetNode();
                INode v1 = (edgeArray[j].getSourceNode().getLayout().getCenter().x <= edgeArray[j].getTargetNode().getLayout().getCenter().x) ? edgeArray[j].getSourceNode() : edgeArray[j].getTargetNode();
                INode v2 = (edgeArray[j].getSourceNode().getLayout().getCenter().x >  edgeArray[j].getTargetNode().getLayout().getCenter().x) ? edgeArray[j].getSourceNode() : edgeArray[j].getTargetNode();

                YPoint p_u1 = new YPoint(u1.getLayout().getCenter().x, u1.getLayout().getCenter().y);
                YPoint p_u2 = new YPoint(u2.getLayout().getCenter().x, u2.getLayout().getCenter().y);
                YPoint p_v1 = new YPoint(v1.getLayout().getCenter().x, v1.getLayout().getCenter().y);
                YPoint p_v2 = new YPoint(v2.getLayout().getCenter().x, v2.getLayout().getCenter().y);

                LineSegment l_e1 = new LineSegment(p_u1, p_u2);
                LineSegment l_e2 = new LineSegment(p_v1, p_v2);

                if (LineSegment.getIntersection(l_e2, l_e1)!=null)
                {
                    numberOfCrossings ++;
                }
                
            }
        }
        return numberOfCrossings;
    }

    public static double calculateCrossingResolution(GraphComponent view)
    {
        double angle = Double.MAX_VALUE;

        IEdge[] edgeArray = new IEdge[view.getGraph().getEdges().size()];
        int k=0;
        for (IEdge e : view.getGraph().getEdges())
        {
            edgeArray[k] = e;
            k++;
        }
        for (int i=0; i<edgeArray.length; i++)
        {
            for (int j=i+1; j<edgeArray.length; j++)
            {
                // e1 = (u1,u2) u1<u2 and e2 = (v1,v2) v1<v2

                INode u1 = (edgeArray[i].getSourceNode().getLayout().getCenter().x <= edgeArray[i].getTargetNode().getLayout().getCenter().x) ? edgeArray[i].getSourceNode() : edgeArray[i].getTargetNode();
                INode u2 = (edgeArray[i].getSourceNode().getLayout().getCenter().x >  edgeArray[i].getTargetNode().getLayout().getCenter().x) ? edgeArray[i].getSourceNode() : edgeArray[i].getTargetNode();
                INode v1 = (edgeArray[j].getSourceNode().getLayout().getCenter().x <= edgeArray[j].getTargetNode().getLayout().getCenter().x) ? edgeArray[j].getSourceNode() : edgeArray[j].getTargetNode();
                INode v2 = (edgeArray[j].getSourceNode().getLayout().getCenter().x >  edgeArray[j].getTargetNode().getLayout().getCenter().x) ? edgeArray[j].getSourceNode() : edgeArray[j].getTargetNode();

                YPoint p_u1 = new YPoint(u1.getLayout().getCenter().x, u1.getLayout().getCenter().y);
                YPoint p_u2 = new YPoint(u2.getLayout().getCenter().x, u2.getLayout().getCenter().y);
                YPoint p_v1 = new YPoint(v1.getLayout().getCenter().x, v1.getLayout().getCenter().y);
                YPoint p_v2 = new YPoint(v2.getLayout().getCenter().x, v2.getLayout().getCenter().y);

                LineSegment l_e1 = new LineSegment(p_u1, p_u2);
                LineSegment l_e2 = new LineSegment(p_v1, p_v2);

                YPoint p = LineSegment.getIntersection(l_e2, l_e1);


                if (p != null)
                {
                    //e1 intersects e2 at p
                    YVector vector_p_u1 = new YVector(p_u1, p);
                    YVector vector_p_u2 = new YVector(p_u2, p);

                    YVector vector_p_v1 = new YVector(p_v1, p);
                    YVector vector_p_v2 = new YVector(p_v2, p);

                    if (YVector.angle(vector_p_u2, vector_p_v2) < Math.PI/2 )
                    {
                        angle = Math.min(angle, YVector.angle(vector_p_u2, vector_p_v2));
                    }
                    else if (YVector.angle(vector_p_v2, vector_p_u2) < Math.PI/2 )
                    {
                        angle = Math.min(angle, YVector.angle(vector_p_v2, vector_p_u2));
                    }
                    else if (YVector.angle(vector_p_v1, vector_p_u2) < Math.PI/2 )
                    {
                        angle = Math.min(angle, YVector.angle(vector_p_v1, vector_p_u2));
                    }
                    else if (YVector.angle(vector_p_u1, vector_p_v2) < Math.PI/2 )
                    {
                        angle = Math.min(angle, YVector.angle(vector_p_u1, vector_p_v2));
                    }
                    else
                    {
                        //javax.swing.JOptionPane.showMessageDialog(null, "unexpected case");
                    }
                }
            }
        }
        return 180*angle/Math.PI;
    }

}
