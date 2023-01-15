/*
 * This file is part of the repicea-simulation library.
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
package repicea.simulation.covariateproviders.plotlevel;

import repicea.util.REpiceaTranslator;
import repicea.util.REpiceaTranslator.TextableEnum;

/**
 * This interface ensures that the stand intance can provide its management type, i.e. either
 * even-aged or uneven-aged.
 * @author Mathieu Fortin - April 2019
 */
public interface ManagementTypeProvider {

	public static enum ManagementType implements TextableEnum {
		UnevenAged("Uneven-aged", "Irr\u00E9gulier"), 
		EvenAged("Even-aged", "R\u00E9gulier");

		ManagementType(String englishText, String frenchText) {
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
	 * This method returns the current management of the stand. 
	 * @return a ManagementType enum
	 */
	public ManagementType getManagementType();

}
