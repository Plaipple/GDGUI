/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author michael
 */

import com.yworks.yfiles.view.GraphComponent;

import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import com.yworks.yfiles.graph.IGraph;
import com.yworks.yfiles.graphml.GraphMLIOHandler;
 
public class Experiment
{
	private GraphComponent view;
    //private y.view.Graph2DView view;

    public Experiment()
    {
    	view = new GraphComponent();
        //this.view = new y.view.Graph2DView();
    }
    
    public void run()
    {
        try
        {
            String csv = "C:/Users/patri_000/Desktop/results.csv";

            java.io.FileWriter fstream = new java.io.FileWriter(csv, false);
            java.io.BufferedWriter out = new java.io.BufferedWriter(fstream);
            out.write("");
            //Close the output stream
            out.close();

            java.text.DecimalFormat df = new java.text.DecimalFormat("###.##");

            java.lang.StringBuffer header = new java.lang.StringBuffer().append("Graph")
                                                                        .append("\t")
                                                                        .append("Nodes")
                                                                        .append("\t")
                                                                        .append("Edges")
                                                                        .append("\t")
                                                                        .append("Max Degree")
                                                                        .append("\t")
                                                                        .append("Eades Angular")
                                                                        .append("\t")
                                                                        .append("Eades Crossing")
                                                                        .append("\t")
                                                                        .append("Eades Iterations")
                                                                        .append("\t")
                                                                        .append("Eades Edge Length")
                                                                        .append("\t")
                                                                        .append("Eades Crossings")
                                                                        .append("\t")
                                                                        .append("Eades Time")
                                                                        .append("\n");

            fstream = new java.io.FileWriter(csv, true);
            out = new java.io.BufferedWriter(fstream);
            out.write(header.toString());
            //Close the output stream
            out.close();

            String directory = "C:/Users/patri_000/Desktop/Graphml/";

            java.io.File dir = new java.io.File(directory);

            // It is also possible to filter the list of returned files.
            java.io.FilenameFilter filter = new java.io.FilenameFilter()
            {
                public boolean accept(java.io.File file, String name)
                {
                    return name.endsWith(".graphml");
                }
            };
            String[] children = dir.list(filter);

            if (children == null)
            {
                // Either dir does not exist or is not a directory
            }
            else
            {
            	GraphMLIOHandler ioh = new GraphMLIOHandler();

                double eadesAngular = 0;
                double eadesCrossing = 0;
                double eadesEdgeLength = 0.0;
                int eadesNoOfCrossings = 0;
                int eadesIterations = 0;
                long eadesTime = 0;

                long startTime = 0;
                long finishTime = 0;

                int maxDegree = 0;

                
                
                for (int i=0; i<children.length; i++)
                {
                    System.out.println("Iteration: " + i);

                    //view.getGraph2D().clear();
                    view.getGraph().clear();
                    //ioh.read(view.getGraph2D(), directory + children[i]);
                    ioh.read(view.getGraph(), directory + children[i]);

                    view.fitContent();
                    view.requestFocus();

                    //maxDegree = util.Utilities.maxDegree(view.getGraph2D());
                    //maxDegree = util.Utilities.maxDegree(view.getGraph());

                   
                    startTime = System.currentTimeMillis();

                    //Then, run Eades' algorithm
                    layout.algo.ForceDirectedAlgorithm eades = new layout.algo.ForceDirectedAlgorithm(view, 1000, 100000, 0.01) {
                        public void calculateVectors() {
                            layout.algo.ForceDirectedFactory.calculateSpringForcesEades(graph, 100, 100, 0.01, map);
                            layout.algo.ForceDirectedFactory.calculateElectricForcesEades(graph, 100000, 0.01, map);
                        }
                    };
                    eades.run();

                    finishTime = System.currentTimeMillis();

                    /*
                    eadesAngular = util.Utilities.calculateAngularResolution(view);
                    eadesCrossing = util.Utilities.calculateCrossingResolution(view);
                    eadesEdgeLength = util.Utilities.calculateAverageEdgeLength(view);
                    eadesNoOfCrossings = util.Utilities.calculateNumberOfCrossings(view);
                    eadesIterations = eades.getMaxNoOfIterations();
                    eadesTime = (finishTime - startTime);
          			*/
                   
                    java.lang.StringBuffer buffer = new java.lang.StringBuffer().append(children[i])
                                                                                .append("\t")
                                                                                //.append(view.getGraph2D().nodeCount())
                                                                                .append(view.getGraph().getNodes().size())
                                                                                .append("\t")
                                                                                //.append(view.getGraph2D().edgeCount())
                                                                                .append(view.getGraph().getEdges().size())
                                                                                .append("\t")
                                                                                .append(maxDegree)
                                                                                .append("\t")
                                                                                .append(df.format(eadesAngular).replace(',', '.'))
                                                                                .append("\t")
                                                                                .append(df.format(eadesCrossing).replace(',', '.'))
                                                                                .append("\t")
                                                                                .append(eadesIterations)
                                                                                .append("\t")
                                                                                .append(df.format(eadesEdgeLength).replace(',', '.'))
                                                                                .append("\t")
                                                                                .append(eadesNoOfCrossings)
                                                                                .append("\t")
                                                                                .append(df.format(eadesTime).replace(',', '.'))
                                                                                .append("\n");

                   
                    fstream = new java.io.FileWriter(csv, true);
                    out = new java.io.BufferedWriter(fstream);
                    out.write(buffer.toString());
                    out.close();
                }
            }
        }
        catch (java.io.IOException exc)
        {
            System.out.println(exc);
        }
    }


    

    public static void main(String[] args)
    {
        Experiment e = new Experiment();
        e.run();
    }
}
