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
public abstract class SimulatedAnnealingAlgorithm  implements Runnable
{
    protected GraphComponent view;
    protected IGraph graph;
    protected int maxNoOfIterations;                //The maximum number of iterations.

    //Graph Listeners.
    // A listener can, e.g., interrupt the iteration process, if it detects converge by setting maxNoOfIterations to -1.
    protected List<AlgorithmListener> algorithmListeners;

    /**
     * Constructor of Objects of type ForceDirectedAlgorithm
     * @param view - an object of type Graph2DView
     * @param maxNoOfIterations - the maximum number of iterations
     */
    public SimulatedAnnealingAlgorithm(GraphComponent view, int maxNoOfIterations)
    {
        this.view = view;
        this.graph = view.getGraph();
        this.maxNoOfIterations = maxNoOfIterations;
        this.algorithmListeners = new ArrayList<AlgorithmListener>();
    }

    /**
     * Abstract method to calculate the vectors.
     * Subclasses must implement this method.
     */
    public abstract void calculatePositions();

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
            this.calculatePositions();
            //this.displayVectors();
        }

        for (int i=0; i<this.maxNoOfIterations; i++)
        {
            this.init();
            this.calculatePositions();
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

