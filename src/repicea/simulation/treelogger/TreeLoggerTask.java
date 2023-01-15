/*
 * This file is part of the repicea-simulation library.
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

import java.util.Collection;

import repicea.app.AbstractGenericTask;

/**
 * The TreeLoggerTask class simply processes the tree logging.
 * @author Mathieu Fortin - January 2013
 */
@SuppressWarnings("serial")
class TreeLoggerTask<Tree extends LoggableTree> extends AbstractGenericTask {

	private TreeLogger<?, Tree> treeLogger;
	
	protected TreeLoggerTask(TreeLogger<?, Tree> treeLogger) {
		this.treeLogger = treeLogger;
	}
	
	@Override
	protected void doThisJob() throws Exception {
		Collection<Tree> loggableTrees = treeLogger.getLoggableTrees();
		if (loggableTrees.size() > 0) {
			double progressFactor = 100d / loggableTrees.size();
			int i = 0;
			setProgress(i);
			for (Tree tree : loggableTrees) {
				if (isCancelled()) {
					break;
				}
				treeLogger.logThisTree(tree);
				setProgress((int) (++i * progressFactor));
			}
		}
	}

}
