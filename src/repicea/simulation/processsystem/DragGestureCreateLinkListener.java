/*
 * This file is part of the repicea-simulation library.
 *
 * Copyright (C) 2009-2014 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation.processsystem;

import java.awt.dnd.DragGestureEvent;
import java.awt.dnd.DragGestureListener;

import repicea.gui.CommonGuiUtility;

public class DragGestureCreateLinkListener implements DragGestureListener {
	
	protected ProcessorButton button;
	
	protected DragGestureCreateLinkListener(ProcessorButton button) {
		this.button = button;
	}
	
	@Override
	public void dragGestureRecognized(DragGestureEvent dge) {
		SystemPanel panel = (SystemPanel) CommonGuiUtility.getParentComponent(button, SystemPanel.class);
		PreProcessorLinkLine futureLink = instantiatePreLink(panel);
		panel.registerLinkBeingCreated(futureLink);
	}
	
	protected PreProcessorLinkLine instantiatePreLink(SystemPanel panel) {
		return new PreProcessorLinkLine(panel, button);
	}
	
}


