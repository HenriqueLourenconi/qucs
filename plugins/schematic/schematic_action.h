/***************************************************************************
    copyright            : (C) 2003 by Michael Margraf
                               2020, 2021 Felix Salfelder
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef QUCS_MOUSEACTION_H
#define QUCS_MOUSEACTION_H

#include "element.h"
#include "action.h"
#include "qt_compat.h"
#include "schematic_scene.h"
#include "element_graphics.h"

class QAction;
class QEvent;
class QMenu;
class QMouseEvent;
class QPainter;
class QUndoCommand;
class Schematic;

namespace qucs {

class MouseActions;

// it is a mode. not an action...
class MouseAction : public QObject, public Action {
public:
	typedef QUndoCommand cmd;
protected:
public:
	explicit MouseAction(QObject*p=nullptr)
		: Action(), _sender(nullptr){(void)p;}
	MouseAction(MouseAction const&) = delete;

public:
	virtual ~MouseAction(){}

public:
	cmd* handle(QEvent*) override;
	cmd* activate(QObject* sender) override;
	cmd* deactivate() override;

// private: TODO
	// TODO: only use POS in those
	virtual cmd* move(QEvent*) { return nullptr; }
	virtual cmd* press(QEvent*) { return nullptr; }
	// virtual cmd* grab(QGraphicsSceneEvent*) { return nullptr; }
	virtual cmd* release(QEvent*) { return nullptr; }
	virtual cmd* dblclk(QEvent*) { return nullptr; }

	virtual cmd* generic(QEvent*) { return nullptr; } // remove
	virtual cmd* enter(QEvent*) {itested(); return nullptr; }
	virtual cmd* leave(QEvent*) {itested(); return nullptr; }

	void uncheck();

protected:
	MouseActions* ctx() const;
	QGraphicsView* view(); // what?
	SchematicScene* scene(); // passed to UndoCommands (check: why?)

protected:
	Doc const* doc() const;
	SchematicScene const* scene() const;

protected:
	void sceneAddItem(ElementGraphics*);
	void sceneRemoveItem(ElementGraphics*);

protected:
	Doc* doc(); // BUG _ctx.
	QList<ElementGraphics*> selectedItems(); // BUG. copies.
	QPointF mapToScene(QPoint const& p) const;
	void updateViewport(); // why?
	void setCursor(QCursor const& c);
   bool isNode(pos_t const&) const; // needed??
   bool isConductor(pos_t const&) const;

protected: // UC
	QList<ElementGraphics*> items(const QPointF &pos,
                                 Qt::ItemSelectionMode mode=Qt::IntersectsItemShape,
                                 Qt::SortOrder order = Qt::DescendingOrder) const;

	Node const* nodeAt(pos_t) const;

private:
//	MouseActions& _ctx; parent?
	QAction* _sender;
}; // MouseAction

} // qucs

#endif