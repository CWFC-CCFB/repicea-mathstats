/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2019 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation;

import repicea.util.REpiceaTranslator;
import repicea.util.REpiceaTranslator.TextableEnum;

/**
 * This interface ensures the stand can provide the scale at which the model applies. It impacts the way
 * the disturbances affect the plots. Typically, the harvesting will be applied to all the plots when the 
 * simulation is run at the stand level. Otherwise, it is asynchronous. 
 * @author Mathieu Fortin - January 2016
 */
public interface ApplicationScaleProvider {
	
	public enum ApplicationScale implements TextableEnum {
		FMU("Regional level", "Echelle r\u00E9gionale"),
		Stand("Stand level", "Echelle du peuplement");
		
		ApplicationScale(String englishText, String frenchText) {
			setText(englishText, frenchText);
		}

		@Override
		public void setText(String englishText, String frenchText) {
			REpiceaTranslator.setString(this, englishText, frenchText);
		}
		
		@Override
		public String toString() {return REpiceaTranslator.getString(this);}
	}

	/**
	 * This method returns the scale at which the model applies. It impacts the disturbances.
	 * @return an ApplicationScale enum
	 */
	public ApplicationScale getApplicationScale();
	
	
}
