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

import java.awt.dnd.DnDConstants;
import java.awt.dnd.DragGestureEvent;
import java.awt.dnd.DragSource;

import javax.swing.Icon;

import repicea.gui.dnd.DragGestureImpl;
import repicea.gui.permissions.REpiceaGUIPermission;

/**
 * This interface ensures a DnDButton in the toolbar provides an appropriate format for the drag and drop.
 * @author Mathieu Fortin - April 2014
 */
@SuppressWarnings("serial")
public abstract class DnDCompatibleButton extends SelectableJButton {

	protected static class InternalDragGestureImpl extends DragGestureImpl<Processor> {
		
		@Override
		protected Processor adaptSourceToTransferable(DragGestureEvent event) {
			Processor obj = ((DnDCompatibleButton) event.getComponent()).createNewProcessor();
			return obj;
		}
	}

	public DnDCompatibleButton(REpiceaGUIPermission permission) {
		super(permission);
		if (permission.isDragGranted()) {
			DragSource ds = new DragSource();
			ds.createDefaultDragGestureRecognizer(this, DnDConstants.ACTION_COPY, new InternalDragGestureImpl());
		}
	}

	/**
	 * This method returns a Processor that is represented by the DnD button in the toolbar.
	 * @return a Processor instance
	 */
	public abstract Processor createNewProcessor();
	
	@Override
	public final void setSelected(boolean bool) {}		// to disable the bold borders

	@Override
	protected Icon getDefaultIcon() {
		return UISetup.Icons.get(ProcessorButton.class.getName());
	}

}
