package layout.algo;

import com.yworks.yfiles.algorithms.YPoint;
import com.yworks.yfiles.algorithms.YVector;
import com.yworks.yfiles.geometry.PointD;
import com.yworks.yfiles.graph.IGraph;
import com.yworks.yfiles.graph.IMapper;
import com.yworks.yfiles.graph.INode;
import com.yworks.yfiles.graph.Mapper;
import com.yworks.yfiles.view.GraphComponent;
import com.yworks.yfiles.view.ICanvasObject;
import com.yworks.yfiles.view.ICanvasObjectDescriptor;
import layout.algo.event.AlgorithmEvent;
import layout.algo.event.AlgorithmListener;
import view.visual.VectorVisual;

import java.awt.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.WeakHashMap;

/**
 * This abstract class is implemented a framework for simulated-annealing algorithms, in which several different
 * criteria can be evaluated through the method calculatePositions(). The implementation is done using the Runnable
 * interface, where AlgorithmListeners are notified during the execution.
 *
 * @author Patrick Laipple
 */
public abstract class SimulatedAnnealingAlgorithm implements Runnable
{
    //Declaration of Variables which are also needed to be saved between the iterations
    protected GraphComponent view;
    protected IGraph graph;
    protected int area;
    protected int temperature;
    protected int index = 0;
    protected double energyOld = 0;
	protected double boundThreshold = 10;  
	protected double bound_top = Double.POSITIVE_INFINITY;
	protected double bound_bottom = Double.NEGATIVE_INFINITY;
	protected double bound_left = Double.POSITIVE_INFINITY;
	protected double bound_right = Double.NEGATIVE_INFINITY;
	protected PointD nodePositions[];
    protected int maxNoOfIterations;                //The maximum number of iterations.

    //Graph Listeners.
    // A listener can, e.g., interrupt the iteration process, if it detects converge by setting maxNoOfIterations to -1.
    protected List<AlgorithmListener> algorithmListeners;

    /**
     * Constructor of Objects of type ForceDirectedAlgorithm
     * @param view - an object of type Graph2DView
     * @param maxNoOfIterations - the maximum number of iterations
     */
    public SimulatedAnnealingAlgorithm(GraphComponent view, int maxNoOfIterations, int area)
    {
        this.view = view;
        this.graph = view.getGraph();
        this.maxNoOfIterations = maxNoOfIterations;
        this.area = area;
        this.algorithmListeners = new ArrayList<AlgorithmListener>();
    }

    /**
     * Abstract method to calculate the Energy Function.
     * Subclasses must implement this method.
     */
    public abstract double calculateEnergyFunction();

    /**
     * Execute the algorithm
     */
    public void run()
    {
        AlgorithmEvent evt = new AlgorithmEvent(this, 0);

        //Notify the listeners that the algorithms is started.
        for (Iterator<AlgorithmListener> it = this.algorithmListeners.iterator(); it.hasNext(); )
        {
            it.next().algorithmStarted(evt);
        }

        // Just for debugging purposes, to display the vectors.
        if (this.maxNoOfIterations == 0)
        {
            this.init();
            //this.calculateEnergyFunction();
            //this.displayVectors();
        }

        
        //Only in the first iteration the bounds are to be set for the area the graph is allowed to use
        //and the array nodePositions gets filled with positions of the nodes in the graph
        nodePositions = new PointD[graph.getNodes().size()];
        int nodeNumber = 0;

        for (INode n : graph.getNodes())
        {
        	n.setTag(nodeNumber);
        	bound_top = Math.min(n.getLayout().getCenter().y, bound_top);
        	bound_bottom = Math.max(n.getLayout().getCenter().y, bound_bottom);
        	bound_left = Math.min(n.getLayout().getCenter().x, bound_left);
        	bound_right = Math.max(n.getLayout().getCenter().x, bound_right);
        	nodePositions[(int)n.getTag()] = new PointD(n.getLayout().getCenter().x, n.getLayout().getCenter().y);
        	nodeNumber ++;
        }

        //The bottom and right bound are set dynamically depending on how many nodes
        //the graph has. The factor 'area' can be set in the panel for the SA Algorithm
        double dynamic_bound_bottom = bound_top + (area * graph.getNodes().size());
        double dynamic_bound_right = bound_left + (area * graph.getNodes().size());
        double graph_center;

        //Don't take the dynamic values if they are greater than the position of the most outside node
        //Because it could happen that nodes of the initial graph are placed outside the dynamic bounds
        //Then these are ignored while running the algorithm and do not change their positions
        if (dynamic_bound_bottom > bound_bottom)
        {
        	//calculate the center of the graph's y-axis. Then add half the value of the 
        	//dynamic space for the bottom and subtract half of it for the top
        	graph_center = bound_top + (bound_bottom - bound_top) / 2;
        	bound_bottom = graph_center + (dynamic_bound_bottom - bound_top) / 2;
        	bound_top = graph_center - (dynamic_bound_bottom - bound_top) / 2;
        }

        if (dynamic_bound_right > bound_right)
        {
        	//same as above except it procedures now for the x-axis
        	graph_center = bound_left + (bound_right - bound_left) / 2;
        	bound_right = graph_center + (dynamic_bound_right - bound_left) / 2;
        	bound_left = graph_center - (dynamic_bound_right - bound_left) / 2;
        }



        
        for (int i=0; i<this.maxNoOfIterations; i++)
        {
            this.init();                      

        	//Adjust the temperature value after each iteration
        	temperature = maxNoOfIterations - index;
        	double newtemperature = ((double)maxNoOfIterations - (double)index) / (double)maxNoOfIterations * 6;
        	double positionRadius;


        	//Calculate the energy function of the initial graph layout
        	energyOld = this.calculateEnergyFunction();

        	//main loop of the Simulated Annealing algorithm
        	for (INode n : graph.getNodes())
        	{   
        		//Energy function of the new graph positioning after one node has been moved
        		double energyNew = 0;   		
        		PointD n_old = nodePositions[(int)n.getTag()];

        		//The radius of the new position is determined dynamically and decreases during the runtime.
        		//(at 1000 iterations it decreases by 0.1%)
        		positionRadius = 100 * ((double)temperature / (double)maxNoOfIterations);
        		//positionRadius = 100 * ((double)(newtemperature * maxNoOfIterations / 6) / (double)maxNoOfIterations);

        		//Creating randomized coordinates for the new position within the distance of positionRadius in any direction
        		int signx = (Math.random() > 0.5) ? -1 : 1;
        		int signy = (Math.random() > 0.5) ? -1 : 1;
        		double newposx = Math.random() * positionRadius * signx;
        		double newposy = Math.random() * positionRadius * signy;
        		PointD n_new = new PointD((n.getLayout().getCenter().x + newposx), (n.getLayout().getCenter().y + newposy));

        		//check if the new position lies within the bounds + some extra space that can be determined manually
        		if(n_new.x > bound_right - boundThreshold || n_new.x < bound_left + boundThreshold || n_new.y > bound_bottom - boundThreshold || n_new.y < bound_top + boundThreshold) continue;

        		//change the position of the node in the array so that the new energy can be calculated for the graph
        		nodePositions[(int)n.getTag()] = n_new;

        		/// Now the same procedure is executed for the new Point
        		energyNew = this.calculateEnergyFunction();
        		
        		//the better positioning leads to saving the coordinates for the node n and updating the actual graph
        		if (energyNew < energyOld)
        		{
        			energyOld = energyNew;
        			graph.setNodeCenter(n, n_new);
        		}
        		//if the new position has a worse energy value the algorithm still might take it
        		//by a probability which is depending on the current temperature and the 
        		//energy difference of the two positions
        		else
        		{
        			/*double ratio = energyOld / energyNew * 0.1;
        			double temperatureRatio = (double) temperature / (double) maxNoOfIterations;
        			double probability = ratio * temperatureRatio;*/
        			double probability = Math.exp(((energyNew - energyOld) * -1) / newtemperature);
        			//System.out.println(probability);
        			if (Math.random() <= probability)
        			{
        				energyOld = energyNew;
        				graph.setNodeCenter(n, n_new);
        			}
        			else
        			{
        				nodePositions[(int)n.getTag()] = n_old;
        			}
        		}    		
        	}

        	index ++;
        	if (temperature <= 1) index = 0;

            
            //this.draw();

            try
            {
                Thread.sleep(1);
            }
            catch (InterruptedException exc)
            {
                //Do nothing...
            }
            this.reset();

            //Notify the listeners that the algorithms changed its status.
            for (Iterator<AlgorithmListener> it = this.algorithmListeners.iterator(); it.hasNext(); )
            {
                evt.currentStatus(Math.round(100*i/this.maxNoOfIterations));
                it.next().algorithmStateChanged(evt);
            }
        }

        //Notify the listeners that the algorithms is finished.
        for (Iterator<AlgorithmListener> it = this.algorithmListeners.iterator(); it.hasNext(); )
        {
            evt.currentStatus(100);
            it.next().algorithmFinished(evt);
        }
    }



	/**
     * Initiates a run.
     */
    protected void init()
    {
    	
    }

    /**
     * Resets the algorithm after each iteration.
     */
    protected void reset()
    {

    }

    /**
     * Adds a new Algorithm Listener.
     * @param algorithmListener - an algorithm listener
     */
    public void addAlgorithmListener(AlgorithmListener algorithmListener)
    {
        this.algorithmListeners.add(algorithmListener);
    }

    /**
     * Remove an Algorithm Listener.
     * @param algorithmListener - an algorithm listener
     */
    public void removeAlgorithmListener(AlgorithmListener algorithmListener)
    {
        this.algorithmListeners.remove(algorithmListener);
    }

    /**
     * Returns the maximum number of iterations
     * @return - the maximum number of iterations
     */
    public int getMaxNoOfIterations() {
        return maxNoOfIterations;
    }

    /**
     * Sets the maximum number of iterations
     * @param maxNoOfIterations - the maximum number of iterations
     */
    public void setMaxNoOfIterations(int maxNoOfIterations) {
        this.maxNoOfIterations = maxNoOfIterations;
    }
}

