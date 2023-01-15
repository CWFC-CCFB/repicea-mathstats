/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2013 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation.covariateproviders.treelevel;

import repicea.util.REpiceaTranslator;
import repicea.util.REpiceaTranslator.TextableEnum;


/**
 * This interface ensures the tree can provide its status.
 * @author Mathieu Fortin - May 2013
 */
public interface TreeStatusProvider {

	
	public static enum StatusClass implements TextableEnum {
		alive("alive", "vivant"), 
		cut("cut", "coup\u00E9"), 
		dead("dead", "mort"),
		windfall("windfall", "chablis");

		StatusClass(String englishText, String frenchText) {
			setText(englishText, frenchText);
		}

		@Override
		public void setText(String englishText, String frenchText) {
			REpiceaTranslator.setString (this, englishText, frenchText);
		}

		@Override
		public String toString() {
			return REpiceaTranslator.getString(this);
		}
	}

	/**
	 * This method sets the status of the tree.
	 * @param statusClass a StatusClass enum
	 */
	public void setStatusClass(StatusClass statusClass);
	
	/**
	 * This method returns the status class of the tree.
	 * @return a StatusClass enum
	 */
	public StatusClass getStatusClass();
	
}
