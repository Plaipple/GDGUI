/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author michael
 */

import com.yworks.yfiles.view.GraphComponent;

import layout.algo.ForceDirectedFactory;
import layout.algo.RandomMovementFactory;
import layout.algo.SimulatedAnnealingFactory;

import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.Duration;

import com.yworks.yfiles.graph.IGraph;
import com.yworks.yfiles.graph.LayoutUtilities;
import com.yworks.yfiles.graphml.GraphMLIOHandler;
import com.yworks.yfiles.layout.GraphTransformer;
import com.yworks.yfiles.layout.OperationType;
import com.yworks.yfiles.layout.organic.OrganicLayout;
 
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
            String txt = "C:/Users/patri_000/Desktop/results.txt";

            java.io.FileWriter fstream = new java.io.FileWriter(txt, false);
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
                                                                        .append("EadesArea")
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
                                                                        .append("Eades Distance")
                                                                        .append("\t")
                                                                        .append("Eades Time")
                                                                        .append("\t")/*
                                                                        .append("\t")
                                                                        .append("Graph")
                                                                        .append("\t")
                                                                        .append("Nodes")
                                                                        .append("\t")
                                                                        .append("Edges")
                                                                        .append("\t")
                                                                        .append("Max Degree")
                                                                        .append("\t")
                                                                        .append("NodeEdgeArea")
                                                                        .append("\t")
                                                                        .append("NodeEdge Angular")
                                                                        .append("\t")
                                                                        .append("NodeEdge Crossing")
                                                                        .append("\t")
                                                                        .append("NodeEdge Iterations")
                                                                        .append("\t")
                                                                        .append("NodeEdge Edge Length")
                                                                        .append("\t")
                                                                        .append("NodeEdge Crossings")
                                                                        .append("\t")
                                                                        .append("NodeEdge Distance")
                                                                        .append("\t")
                                                                        .append("NodeEdge Time")
                                                                        .append("\n")*/
                                                                        .append("\t")
                                                                        .append("Graph")
                                                                        .append("\t")
                                                                        .append("Nodes")
                                                                        .append("\t")
                                                                        .append("Edges")
                                                                        .append("\t")
                                                                        .append("Max Degree")
                                                                        .append("\t")
                                                                        .append("RM Area")
                                                                        .append("\t")
                                                                        .append("RM Angular")
                                                                        .append("\t")
                                                                        .append("RM Crossing")
                                                                        .append("\t")
                                                                        .append("RM Iterations")
                                                                        .append("\t")
                                                                        .append("RM Edge Length")
                                                                        .append("\t")
                                                                        .append("RM Crossings")
                                                                        .append("\t")
                                                                        .append("RM Distance")
                                                                        .append("\t")
                                                                        .append("RM Time")
                                                                        .append("\n");

            fstream = new java.io.FileWriter(txt, true);
            out = new java.io.BufferedWriter(fstream);
            out.write(header.toString());
            //Close the output stream
            out.close();

            String inputDirectory = "C:/Users/patri_000/Desktop/Graphml/";
            String outputDirectoryEades = "C:/Users/patri_000/Desktop/OutputEades/";
            String outputDirectoryNodeEdge = "C:/Users/patri_000/Desktop/OutputSA/";
            String outputDirectoryRandomMovement = "C:/Users/patri_000/Desktop/OutputRandMov/";

            java.io.File dir = new java.io.File(inputDirectory);

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
                double eadesDistance = 0.0;
                double eadesArea = 0.0;
                int eadesNoOfCrossings = 0;
                int eadesIterations = 0;
                long eadesTime = 0;
                
                double randMovAngular = 0;
                double randMovCrossing = 0;
                double randMovEdgeLength = 0.0;
                double randMovDistance = 0.0;
                double randMovArea = 0.0;
                int randMovNoOfCrossings = 0;
                int randMovIterations = 0;
                long randMovTime = 0;

                long startTime = 0;
                long finishTime = 0;

                int maxDegree = 0;

                
                
                for (int i=0; i<children.length; i++)
                {
                    System.out.println("Iteration: " + i);

                    //view.getGraph2D().clear();
                    view.getGraph().clear();
                    //ioh.read(view.getGraph2D(), directory + children[i]);
                    ioh.read(view.getGraph(), inputDirectory + children[i]);

                    view.fitContent();
                    view.requestFocus();

                    maxDegree = util.Utilities.maxDegree(view);

                   
                    startTime = System.currentTimeMillis();

                    //Then, run Eades' algorithm
                    /*
                    layout.algo.ForceDirectedAlgorithm eades = new layout.algo.ForceDirectedAlgorithm(view, 1000) {
                        public void calculateVectors() {
                            layout.algo.ForceDirectedFactory.calculateSpringForcesEades(graph, 150, 150, 0.01, map);
                            layout.algo.ForceDirectedFactory.calculateElectricForcesEades(graph, 100000, 0.01, map);
                            layout.algo.ForceDirectedFactory.calculateElectricForcesCrossingResolution(graph, 100000, 0.01, map);
                            layout.algo.ForceDirectedFactory.calculateElectricForcesAngularResolution(graph, 100000, 0.01, map);
                            //layout.algo.ForceDirectedFactory.calculateElectricForcesNodeEdge(graph, 100000, 0.01, map);
                        }
                    };
                    
                    LayoutUtilities.applyLayout(view.getGraph(), new OrganicLayout());
                                       
                    eades.run();

                    finishTime = System.currentTimeMillis();
                    
                    ioh.write(view.getGraph(), outputDirectoryEades + children[i]);
                    
                    eadesAngular = util.Utilities.calculateAngularResolution(view);
                    eadesCrossing = util.Utilities.calculateCrossingResolution(view);
                    eadesEdgeLength = util.Utilities.calculateAverageEdgeLength(view);
                    eadesNoOfCrossings = util.Utilities.calculateNumberOfCrossings(view);
                    eadesDistance = util.Utilities.calculateshortestNodeEdgeDistance(view);
                    eadesArea = util.Utilities.calculateUsedArea(view);
                    eadesIterations = eades.getMaxNoOfIterations();
                    eadesTime = (finishTime - startTime);*/
                    
                    view.getGraph().clear();
                    ioh.read(view.getGraph(), inputDirectory + children[i]);
                    view.fitContent();
                    view.requestFocus();
                    startTime = System.currentTimeMillis();
                    //Then, run the RandomMovement algorithm
                    layout.algo.RandomMovementAlgorithm randMov= new layout.algo.RandomMovementAlgorithm(view, 1000) {
                        public void calculatePositions() {
                        	RandomMovementFactory.randomMovement(graph, 8, 10, 250, true, true, 10, 10, false, false, false, true, false, false, maxNoOfIterations); 
                        }
                    };
          			
                    LayoutUtilities.applyLayout(view.getGraph(), new OrganicLayout());
                    
                    GraphTransformer gt = new GraphTransformer();
                    gt.setOperation(OperationType.SCALE);
                    gt.setScaleFactor(Double.valueOf(3));
                    LayoutUtilities.morphLayout(view, gt, Duration.ofSeconds(0), null);
                    
                    randMov.run();

                    finishTime = System.currentTimeMillis();
                    
                    ioh.write(view.getGraph(), outputDirectoryRandomMovement + children[i]);
                    
                    randMovAngular = util.Utilities.calculateAngularResolution(view);
                    randMovCrossing = util.Utilities.calculateCrossingResolution(view);
                    randMovEdgeLength = util.Utilities.calculateAverageEdgeLength(view);
                    randMovNoOfCrossings = util.Utilities.calculateNumberOfCrossings(view);
                    randMovDistance = util.Utilities.calculateshortestNodeEdgeDistance(view);
                    randMovArea = util.Utilities.calculateUsedArea(view);
                    randMovIterations = randMov.getMaxNoOfIterations();
                    randMovTime = (finishTime - startTime);
                   
                    java.lang.StringBuffer buffer = new java.lang.StringBuffer().append(children[i])
                                                                                .append("\t")
                                                                                .append(view.getGraph().getNodes().size())
                                                                                .append("\t")
                                                                                .append(view.getGraph().getEdges().size())
                                                                                .append("\t")
                                                                                .append(maxDegree)
                                                                                .append("\t")
                                                                                .append(df.format(eadesArea))//.replace(',', '.'))
                                                                                .append("\t")
                                                                                .append(df.format(eadesAngular))//.replace(',', '.'))
                                                                                .append("\t")
                                                                                .append(df.format(eadesCrossing))//.replace(',', '.'))
                                                                                .append("\t")
                                                                                .append(eadesIterations)
                                                                                .append("\t")
                                                                                .append(df.format(eadesEdgeLength))//.replace(',', '.'))
                                                                                .append("\t")
                                                                                .append(eadesNoOfCrossings)
                                                                                .append("\t")
                                                                                .append(df.format(eadesDistance))//.replace(',', '.'))
                                                                                .append("\t")
                                                                                .append(df.format(eadesTime))//.replace(',', '.'))
                                                                                .append("\t")
                                                                                .append("\t")
                                                                                .append(children[i])
                                                                                .append("\t")
                                                                                .append(view.getGraph().getNodes().size())
                                                                                .append("\t")
                                                                                .append(view.getGraph().getEdges().size())
                                                                                .append("\t")
                                                                                .append(maxDegree)
                                                                                .append("\t")
                                                                                .append(df.format(randMovArea))//.replace(',', '.'))
                                                                                .append("\t")
                                                                                .append(df.format(randMovAngular))//.replace(',', '.'))
                                                                                .append("\t")
                                                                                .append(df.format(randMovCrossing))//.replace(',', '.'))
                                                                                .append("\t")
                                                                                .append(randMovIterations)
                                                                                .append("\t")
                                                                                .append(df.format(randMovEdgeLength))//.replace(',', '.'))
                                                                                .append("\t")
                                                                                .append(randMovNoOfCrossings)
                                                                                .append("\t")
                                                                                .append(df.format(randMovDistance))//.replace(',', '.'))
                                                                                .append("\t")
                                                                                .append(df.format(randMovTime))//.replace(',', '.'))
                                                                                .append("\n");

                   
                    fstream = new java.io.FileWriter(txt, true);
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
