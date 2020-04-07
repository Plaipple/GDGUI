/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package util;

import com.yworks.yfiles.algorithms.YVector;
import com.yworks.yfiles.graph.IEdge;
import com.yworks.yfiles.graph.IGraph;
import com.yworks.yfiles.graph.INode;

/**
 *
 * @author fouli
 */
public class CyclicEdgeComparator implements java.util.Comparator<IEdge>
{
    private IGraph graph;

    public CyclicEdgeComparator(IGraph graph)
    {
        this.graph = graph;
    }

    /**
     * Compares two edges that must share a common end point.
     */
    public int compare(IEdge e1, IEdge e2)
    {
        INode c;
        INode u;
        INode v;

        if (e1.getSourceNode() == e2.getSourceNode())
        {
            c = e1.getSourceNode();
            u = e1.getTargetNode();
            v = e2.getTargetNode();
        }
        else if (e1.getSourceNode() == e2.getTargetNode())
        {
            c = e1.getSourceNode();
            u = e1.getTargetNode();
            v = e2.getSourceNode();
        }
        else if (e1.getTargetNode() == e2.getSourceNode())
        {
            c = e1.getTargetNode();
            u = e1.getSourceNode();
            v = e2.getTargetNode();
        }
        else if (e1.getTargetNode() == e2.getTargetNode())
        {
            c = e1.getTargetNode();
            u = e1.getSourceNode();
            v = e2.getSourceNode();
        }
        else
        {
            return -1;
        }

        YVector cVector = new YVector(c.getLayout().getCenter().x+1, c.getLayout().getCenter().y, c.getLayout().getCenter().x, c.getLayout().getCenter().y);
        YVector uVector = new YVector(u.getLayout().getCenter().x, u.getLayout().getCenter().y, c.getLayout().getCenter().x, c.getLayout().getCenter().y);
        YVector vVector = new YVector(v.getLayout().getCenter().x, v.getLayout().getCenter().y, c.getLayout().getCenter().x, c.getLayout().getCenter().y);

        double tu = YVector.angle(uVector, cVector);
        double tv = YVector.angle(vVector, cVector);

        if (tu == tv)
        {
            return 0;
        }
        else if (tu > tv)
        {
            if (e1.getSourceNode() == u || e1.getTargetNode() == u)
            {
                return 1;
            }
            else
            {
                return -1;
            }
        }
        else
        {
            if (e1.getSourceNode() == u || e1.getTargetNode() == u)
            {
                return -1;
            }
            else
            {
                return 1;
            }
        }
    }
}
