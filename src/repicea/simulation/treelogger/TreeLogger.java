/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2012 Mathieu Fortin for Rouge-Epicea
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * This library is distributed with the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * Please see the license at http://www.gnu.org/copyleft/lesser.html.
 */
package repicea.simulation.treelogger;

import java.awt.Window;
import java.beans.PropertyChangeListener;
import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import repicea.app.GenericTask;
import repicea.gui.genericwindows.REpiceaProgressBarDialog;

/**
 * The TreeLogger abstract class is the class from which a the tree loggers must inherit. Among others,
 * this class contains <br>
 * <br>
 * <li> a Collection of LoggableTree instances that defines the trees to be processed;</li>
 * <li> a SetMap that contains the LoggableTree as key and the resulting WoodPiece instances as values; </li> 
 * <li> a TreeLoggerParameter instance which sets the parameter of the logger.</li><br>
 * <br>
 * The class can be used either in GUI or in script mode this way:<br>
 * <br>
 * {@code MyTreeLogger treeLogger = new MyTreeLogger();} </br>
 * {@code treeLogger.init(MyCollectionOfLoggableTreeInstances);} </br>
 * <br>
 * {@code treeLogger.setTreeLoggerParameters(MyTreeLoggerParameters);		// in script mode} </br>
 * or<br>
 * {@code treeLogger.setTreeLoggerParameters();		// in Gui mode, a dialog will come out} </br>
 * <br>
 * and to process the trees...<br>
 * <br>
 * {@code treeLogger.run(); } <br>
 * <br>
 * This class implements the logging process in a TreeLoggerTask instance which extends the SwingWorker class and can fire events 
 * to eventual listeners. The return type is simply object. If the thread ended
 * correctly, the get() method returns TreeLogger.CORRECTLY_TERMINATED. Otherwise, it returns an exception.<br>
 * <br>
 * @param <Parameter> the class that defines the parameters of this tree logger. It must be a TreeLoggerParameters instance.
 * @author Mathieu Fortin - April 2010
 */
public abstract class TreeLogger<Parameter extends TreeLoggerParameters<? extends LogCategory>, Tree extends LoggableTree> implements GenericTask {
	
	
	@Deprecated
	protected TreeLoggerWrapper wrapper;
	private boolean saveMemory;
	private final Collection<Tree> trees;
	private Map<LoggableTree, Collection<WoodPiece>> woodPieces;
	protected Parameter params;
	protected TreeLoggerTask<?> loggerTask;
	protected List<PropertyChangeListener> listeners;
	
	protected Window owner = null;
	protected boolean progressBarEnabled;
	
	/**
	 * General construtor for all AbstractTreeLogger-derived classes.
	 */
	protected TreeLogger() {
		woodPieces = new HashMap<LoggableTree, Collection<WoodPiece>>();
		listeners = new ArrayList<PropertyChangeListener>();
		setSaveMemoryEnabled(true);		// default value
		trees = new ArrayList<Tree>();
	}

	/**
	 * This method adds the listener instance if it is not already contained in the listeners list.
	 * @param listener a PropertyChangeListener
	 */
	public void addTreeLoggerListener(PropertyChangeListener listener) {
		if (!listeners.contains(listener)) {
			listeners.add(listener);
		}
	}
	
	/**
	 * This method removes the instance listener from the listeners.
	 * @param listener a PropertyChangeListener instance
	 */
	public void removeTreeLoggerListener(PropertyChangeListener listener) {
		while (listeners.contains(listener)) {
			listeners.remove(listener);
		}
	}
	
	
	
	/**
	 * This method initializes the TreeLogger instance by setting the collection of LoggableTree objects 
	 * on which the logger has to work. If the number of LoggableTree objects exceeds the value of 
	 * TreeLogger.NB_TREES_BEFORE_MEMORY_OVERFLOW, the logger adopts a safe memory strategy, i.e., it calls the
	 * method setSaveMemoryEnabled(true). 
	 * @param newtrees a Collection of LoggableTree objects to be processed
	 * @throws InvalidParameterException if the collection trees is null or empty
	 */
	public void init(Collection<?> newtrees) {
		trees.clear();
		trees.addAll(getEligibleTreesFromACollection(newtrees));
	}

	
	/**
	 * This method returns the collection of LoggableTree instances that the logger has to work on. 
	 * @return a Collection of LoggableTree-derived instances
	 */
	public Collection<Tree> getLoggableTrees() {return trees;}
	
	
	/**
	 * This method enables a progress bar related to window owner.
	 * @param bool true to enable the progress bar or false otherwise (which is the default value)
	 * @param owner a Window instance (can be null)
	 */
	public void setProgressBarEnabled(boolean bool, Window owner) {
		this.progressBarEnabled = bool;
		this.owner = owner;
	}
	
	/**
	 * This method adds a WoodPiece instance in the set of WoodPiece associated to the LoggableTree tree.
	 * @param tree a LoggableTree
	 * @param pieceToAdd a WoodPiece to add
	 */
	protected void addWoodPiece(Tree tree, WoodPiece pieceToAdd) {
		Collection<WoodPiece> oColl = woodPieces.get(tree);
		if (oColl == null) {
			oColl = new ArrayList<WoodPiece>();
			woodPieces.put(tree, oColl);
		}
		oColl.add(pieceToAdd);
	}
	
	/**
	 * This method process the trees into logs. It has to be defined in the derived class.
	 */
	protected abstract void logThisTree(Tree tree);

	@Override
	public final void run() {
		woodPieces.clear();
		logTrees();
	}

	/**
	 * This method instantiate a TreeLoggerTask and run this task. Instructions prior to this can be implemented
	 * in the priorToRunning() method while instructions after the completion of this task can be
	 * implemented in the posteriorToRunning() method.
	 */
	@SuppressWarnings("deprecation")
	protected void logTrees() {
		priorToRunning();
		loggerTask = new TreeLoggerTask<Tree>(this);
		for (PropertyChangeListener listener : listeners) {
			loggerTask.addPropertyChangeListener(listener);
		}
		if (progressBarEnabled) {
			new REpiceaProgressBarDialog(owner, "title", "message", loggerTask, false);
		} else {
			loggerTask.run();
		}
		for (PropertyChangeListener listener : listeners) {
			loggerTask.removePropertyChangeListener(listener);
		}
		if (isCorrectlyTerminated()) {
			posteriorToRunning();
		} 
	}
	
	/**
	 * This method is run just before launching the logger task.
	 */
	protected void priorToRunning() {}
	

	/**
	 * This method is run just after launching the logger task.
	 */
	protected void posteriorToRunning() {}

	
	/**
	 * This method returns the SetMap object that contains the trees in keys and their 
	 * associated wood pieces once they are processed.
	 * @return a Map<LoggableTree, Collection<WoodPiece>> instance
	 */
	public Map<LoggableTree, Collection<WoodPiece>> getWoodPieces() {return woodPieces;}

	
	/**
	 * This method returns the parameters of the TreeLogger object.
	 * @return a Parameter-defined instance (see header of this class)
	 */
	public Parameter getTreeLoggerParameters() {return params;}

	/**
	 * This method enables the save memory mode. Might have no effect if the tree logger does not
	 * implement any different method. By default this mode is enabled.
	 * @param saveMemory true to enable this mode
	 */
	public void setSaveMemoryEnabled(boolean saveMemory) {this.saveMemory = saveMemory;}
	
	
	/**
	 * This method checks if the logger is operating in safe memory mode.
	 * @return a boolean
	 */
	protected boolean isSaveMemoryEnabled() {return this.saveMemory;}
	

	
	/**
	 * This method sets the parameters in script mode.
	 * @param params a TreeLoggerParameters instance
	 */
	public void setTreeLoggerParameters(Parameter params) {
		if (params != null) {
			if (!params.isCorrect()) {
				throw new InvalidParameterException("The treelogger parameters are incorrect. Please check.");
			}
			this.params = params;
			params.treeLogger = this;
		}
	}
	
	/**
	 * This method is called in GUI mode.
	 */
	public abstract void setTreeLoggerParameters();

	/**
	 * This method returns a set of default parameters to make sure the logger can work. 
	 * @return a set of Parameters-derived instance
	 */
	public abstract Parameter createDefaultTreeLoggerParameters();
	

//	/**
//	 * This static method makes it possible to extract the LoggableTree objects from
//	 * a collection of objects.
//	 * @param trees a Collection of Object instances
//	 * @return a Collection of LoggableTree instances
//	 */
//	public static Collection<? extends LoggableTree> findLoggableTrees(Collection<?> trees) {
//      	Collection<LoggableTree> loggableTrees = new ArrayList<LoggableTree>();
//      	for (Object t : trees) {
//      		if (t instanceof LoggableTree) {
//      			loggableTrees.add((LoggableTree) t);
//      		}
//      	}
//      	return loggableTrees;
//	}
	
	@Override
	public String getName() {
		return "TreeLoggerTask";
	}
	
	@Override
	public void cancel() {
		if (loggerTask != null) {
			loggerTask.cancel();
		}
	}

	@Override
	public boolean isCancelled() {
		if (loggerTask != null) {
			return loggerTask.isCancelled();
		} else {
			return false;
		}
	}

	@Override
	public boolean isCorrectlyTerminated() {
		if (loggerTask != null) {
			return loggerTask.isCorrectlyTerminated();
		} else {
			return false;
		}
	}

	/**
	 * This method extracts a collection of TreePetroLoggable objects from a collection of LoggableTree instances.
	 * @param trees a Collection of LoggableTree-derived instances
	 * @return a Collection of PetroLoggableTree instances
	 */
	protected Collection<Tree> getEligibleTreesFromACollection(Collection<?> trees) {
		Collection<Tree> eligibleTrees = new ArrayList<Tree>();
		if (trees != null) {
			for (Object t : trees) {
				if (t instanceof LoggableTree) {
					Tree tree = getEligible((LoggableTree) t);
					if (tree != null) {
						eligibleTrees.add(tree);
					}
				}
			}
		}
		return eligibleTrees;
	}


	@Override
	public Exception getFailureReason() {
		if (loggerTask != null) {
			return loggerTask.getFailureReason();
		} else {
			return null;
		}
	}

	/**
	 * This method determines whether or not the tree is eligible for this tree logger. Typically, some criteria based on dbh or height can be 
	 * specified here.
	 * @param t a LoggableTree instance
	 * @return the appropriate Tree instance if eligible or null otherwise
	 */
	public abstract Tree getEligible(LoggableTree t);
	
	/**
	 * This method makes it possible to determine whether or not this TreeLogger
	 * class is compatible with the reference object
	 * @param check a TreeLoggerCompatibilityCheck instance
	 * @return a boolean
	 */
	public abstract boolean isCompatibleWith(TreeLoggerCompatibilityCheck check);
	
	@Override
	public boolean isVerbose() {return false;}

		
}
