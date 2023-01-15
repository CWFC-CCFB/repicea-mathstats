/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2015 Mathieu Fortin for Rouge-Epicea
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
package repicea.math.optimizer;

/**
 * This interface should be implemented by AbstractMathematicalFunction-derived classes in which
 * the different steps of the optimization has to be taken into account.
 * @author Mathieu Fortin - November 2015
 */
public interface OptimizerListener {
	
	public static String optimizationStarted = "OptimizerStarted";
	public static String optimizationEnded = "OptimizerEnded";
	
	public void optimizerDidThis(String actionString);
	
}
