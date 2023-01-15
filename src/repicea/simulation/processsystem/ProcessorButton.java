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

import java.awt.Container;
import java.awt.Dimension;
import java.awt.Image;
import java.awt.Point;
import java.awt.Window;
import java.awt.dnd.DnDConstants;
import java.awt.dnd.DragGestureRecognizer;
import java.awt.dnd.DragSource;
import java.awt.event.MouseEvent;

import javax.swing.Icon;
import javax.swing.ImageIcon;

import repicea.gui.CommonGuiUtility;
import repicea.gui.REpiceaAWTProperty;
import repicea.gui.REpiceaShowableUIWithParent;
import repicea.gui.dnd.DragGestureMoveComponentHandler;
import repicea.gui.permissions.REpiceaGUIPermission;
import repicea.gui.permissions.REpiceaGUIPermissionProvider;
import repicea.simulation.processsystem.UISetup.BasicMode;

@SuppressWarnings("serial")
public class ProcessorButton extends SelectableJButton implements AnchorProvider, REpiceaShowableUIWithParent, REpiceaGUIPermissionProvider {

	
	protected static class DragGestureButtonMoveHandler extends DragGestureMoveComponentHandler<ProcessorButton> {
		
		private Point currentCorner;
		private boolean isDragging = false;
		
		protected DragGestureButtonMoveHandler(ProcessorButton component) {
			super(component);
		}

		@Override
		public void mouseDragged(MouseEvent arg0) {
			isDragging = true;
			Point cursorPoint = arg0.getPoint();
			currentCorner = component.getLocation();
			int xCoord = cursorPoint.x - offsetPoint.x + currentCorner.x;
			int yCoord = cursorPoint.y - offsetPoint.y + currentCorner.y;
			Container container = component.getParent();
			if (container != null) {
				if (xCoord + component.getSize().width > container.getSize().width) {
					xCoord = container.getSize().width - component.getSize().width;
				}
				if (xCoord < 0) {
					xCoord = 0;
				}
				if (yCoord + component.getSize().height > container.getSize().height) {
					yCoord = container.getSize().height - component.getSize().height;
				}
				if (yCoord < 0) {
					yCoord = 0;
				}
			}
			Point newLocation = new Point(xCoord, yCoord);
			component.setLocation(newLocation);	
			component.validate();
		}
		
		
		@Override
		public void mouseReleased(MouseEvent evt) {
			if (isEnabled()) {
				super.mouseReleased(evt);
				SystemManagerDialog dlg = (SystemManagerDialog) CommonGuiUtility.getParentComponent(component, SystemManagerDialog.class);
				if (dlg != null && isDragging) {
					dlg.firePropertyChange(REpiceaAWTProperty.ActionPerformed, offsetPoint, component.getLocation());
					isDragging = false;
				}
			}
		}
	}

	
	private final Processor process;
	protected final DragGestureMoveComponentHandler<ProcessorButton> buttonMoveRecognizer;
	protected final DragGestureRecognizer createLinkRecognizer;
	protected final SystemLabel label;
	private final SystemPanel panel;
	private boolean hasChanged;
	protected transient ProcessorInternalDialog guiInterface;
	
	protected ProcessorButton(SystemPanel panel, Processor process, REpiceaGUIPermission permission) {
		super(permission);
		this.panel = panel;
		this.process = process;
		this.label = new SystemLabel();
		addMouseListener(new SystemComponentMouseAdapter(this));
		DragSource ds = new DragSource();
		buttonMoveRecognizer = new DragGestureButtonMoveHandler(this);
		buttonMoveRecognizer.setEnabled(false);
		createLinkRecognizer = ds.createDefaultDragGestureRecognizer(this, DnDConstants.ACTION_COPY, new DragGestureCreateLinkListener(this));
		createLinkRecognizer.setComponent(null);
		setLabel();
	}

	protected ProcessorButton(SystemPanel panel, Processor process) {
		this(panel, process, panel.getListManager().getGUIPermission());
	}
	
	protected void setLabel() {
		String processName = process.getName();
		if (processName.length() > 30) {
			processName = processName.substring(0, 30) + "...";
		}
		label.setText(processName);
	}
	

	@Override
	protected void finalize() {
		panel.remove(this);
		panel.remove(label);
	}
	
	@Override
	public Icon getIcon() {
		if (originalIcon == null) {
			return originalIcon;
		} else{
			Image newImage = originalIcon.getImage().getScaledInstance(SystemLayout.convertToRelative(originalIcon.getIconWidth()),
					SystemLayout.convertToRelative(originalIcon.getIconHeight()),
					Image.SCALE_DEFAULT);
			return new ImageIcon(newImage);
		}
	}


	/**
	 * This method returns the Processor instance that owns this button.
	 * @return a Processor instance
	 */
	public Processor getOwner() {return process;}
		
	@Override
	public Point getRightAnchor() {
		int x = getLocation().x + getSize().width;
		int y = getLocation().y + getSize().height / 2;
		return new Point(x,y);
	}

	@Override
	public Point getLeftAnchor() {
		int x = getLocation().x;
		int y = getLocation().y + getSize().height / 2;
		return new Point(x,y);
	}

	@Override
	public Point getTopAnchor() {
		int x = getLocation().x + getSize().width / 2;
		int y = getLocation().y;
		return new Point(x,y);
	}

	@Override
	public Point getBottomAnchor() {
		int x = getLocation().x + getSize().width / 2;
		int y = getLocation().y + getSize().height;
		return new Point(x,y);
	}

	protected Point retrieveFormerLocation() {
		return process.getOriginalLocation();
	}
	
	@Override
	public void setLocation(int x, int y) {
		super.setLocation(x, y);
		getOwner().setOriginalLocation(new Point(x,y));
		label.setLocation(x, y + (int) getSize().getHeight());
	}

	@Override
	public void setLocation(Point point) {setLocation(point.x, point.y);}
	
	@Override
	public Dimension getPreferredSize() {
		return SystemLayout.convertOriginalToRelative(super.getPreferredSize());
	}
	
	@SuppressWarnings("rawtypes")
	protected void setDragMode(Enum mode) {
		if (getGUIPermission().isDragGranted()) {
			buttonMoveRecognizer.setEnabled(false);
			createLinkRecognizer.setComponent(null);
			if (mode == BasicMode.MoveProcessor) {
				buttonMoveRecognizer.setEnabled(true);
			} else if (mode == BasicMode.CreateLink) {
				if (!getOwner().isTerminalProcessor()) {
					createLinkRecognizer.setComponent(this);
				}
			}
		}
	}

	@Override
	public ProcessorInternalDialog getUI(Container parent) {
		if (guiInterface == null) {
			guiInterface = new ProcessorInternalDialog((Window) parent, this);
		}
		return guiInterface;
	}

	@Override
	public void showUI(Window window) {
		ProcessorInternalDialog dlg = getUI(window);
		if (!dlg.isVisible()) {
			hasChanged = false;
		}
		dlg.setVisible(true);
		if (hasChanged) {
			SystemManagerDialog systemDlg = (SystemManagerDialog) CommonGuiUtility.getParentComponent(this, SystemManagerDialog.class);
			if (systemDlg != null) {
				systemDlg.firePropertyChange(REpiceaAWTProperty.ActionPerformed, this, null);
			}
			hasChanged = false;
		}
	}

	@Override
	public String toString() {
		if (getOwner() != null) {
			return getOwner().getName();
		} else {
			return "";
		}
	}


	@Override
	protected Icon getDefaultIcon() {return UISetup.Icons.get(ProcessorButton.class.getName());}

	protected void setChanged(boolean hasChanged) {this.hasChanged = hasChanged;}
}
