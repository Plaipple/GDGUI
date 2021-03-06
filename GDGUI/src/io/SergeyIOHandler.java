package io;

import com.yworks.yfiles.algorithms.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Created by Michael on 25/6/2016.
 */
public class SergeyIOHandler extends GraphIOHandler {
    /**
     * Reads the contents of the given file and returns a book embedding
     * @param inputFileName the input file
     * @return the graph
     * @throws java.io.IOException
     */
    public static Graph read(String inputFileName) throws IOException
    {
        String line;
        StringBuffer content = new StringBuffer();
        BufferedReader reader = new BufferedReader(new FileReader(inputFileName));
        while((line = reader.readLine()) != null)
        {
            content.append(line).append("\n");
        }
        reader.close();

        Graph g = new Graph();
        Node[] nodes = new Node[0];
        int[][] edges = new int[0][0];

        String[] lines = content.toString().split("\\n+");
        for (int l =0; l < lines.length; l++) {
            String[] source = lines[l].split("\\s+");
            if (l == 0) {
                nodes = new Node[Integer.parseInt(source[0])];
                edges = new int[Integer.parseInt(source[0])][Integer.parseInt(source[0])];
                for (int i = 0; i < nodes.length; i++) {
                    nodes[i] = g.createNode();
                }
            }
            else
            {

                String[] target = source[1].split(",");
                for (int i = 0; i < target.length; i++)
                {
                    //Symmetry
                    edges[Integer.parseInt(source[0])][Integer.parseInt(target[i])] = 1;
                    edges[Integer.parseInt(target[i])][Integer.parseInt(source[0])] = 1;
                }
            }
        }
        for (int i=0; i<edges.length; i++)
        {
            for (int j=edges.length-1; j>i; j--)
            {
                if (edges[i][j] == 1)
                {
                    g.createEdge(nodes[i], nodes[j]);
                }
            }
        }
        return g;
    }

    /**
     * Writes the contents of the given graph to a stream
     * @param graph the graph
     * @param outputFileName the output file
     * @throws java.io.IOException
     */
    public static void write(Graph graph, String outputFileName) throws IOException
    {
        StringBuffer output = new StringBuffer();
        output.append(graph.nodeCount() + " " + graph.edgeCount() + "\n");
        INodeMap map = graph.createNodeMap();
        int k = 0;
        for (INodeCursor nc = graph.getNodeCursor(); nc.ok(); nc.next())
        {
            map.setInt(nc.node(), k);
            k++;
        }
        boolean[][] edges = new boolean[graph.nodeCount()][graph.nodeCount()];
        for (IEdgeCursor ec = graph.getEdgeCursor(); ec.ok(); ec.next())
        {
            edges[map.getInt(ec.edge().source())][map.getInt(ec.edge().target())] = true;
            edges[map.getInt(ec.edge().target())][map.getInt(ec.edge().source())] = true;
        }

        for (int i=0; i<edges.length; i++)
        {
            output.append(i + " ");
            boolean edgeExists = false;
            for (int j=0; j<edges.length; j++)
            {
                if (edges[i][j]) {
                    output.append(j + "," );
                    edgeExists = true;
                }
            }
            if (edgeExists) {
                output.setLength(output.length() - 1);
            }
            output.append("\n");
        }
        PrintWriter out = new PrintWriter(outputFileName);
        out.println(output);
        out.close();
    }
}
