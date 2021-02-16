/***************************************************************************
    copyright            : (C) 2003, 2004 by Michael Margraf
                               2019 Felix Salfelder / QUCS team
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <QFileInfo>
#include <QAction>

#include "action.h"
#include "qucsdoc.h"
#include "doc_actions.h"
#include "simmessage.h"
#include "qucs_app.h"
#include <QUndoStack>
#include "qio.h" // tmp?
#include "simulator.h"
#include "sckt_base.h"
#include "qt_compat.h"

void QucsDoc::setDocFile(std::string const& filename)
{ untested();
	_name = QString_(filename);
  QFileInfo Info(_name);
  if(!_name.isEmpty()) { untested();
    _name = Info.absoluteFilePath();

	 // properties.set("filename", _name); ??
    QString base = Info.completeBaseName();
    QString ext = Info.suffix();

    if(ext == "m" || ext == "oct")
      SimTime = "1";

    DataSet = base + ".dat";       // name of the default dataset
    Script = base + ".m";          // name of the default script
    if(ext != "dpl")
      DataDisplay = base + ".dpl"; // name of default data display
    else { untested();
      DataDisplay = base + ".sch"; // name of default schematic
      GridOn = false;              // data display without grid (per default)
    }
  }else{ untested();
  }
}

QucsDoc::QucsDoc() : Widget() // (QucsApp* App_, const QString& Name_, QWidget* o)
//   : _app(App_),
//	  _owner(o)
{ untested();
  GridOn = true;
  _name = "unnamed"; // name? filename?

  // setDocName(...)?

  SimOpenDpl = true;
  SimRunScript = false;

  DocChanged = false;
  showBias = -1;  // don't show DC bias (currently for "Schematic" only)
  // Scale = 1.0;
}

QucsDoc::~QucsDoc()
{ untested();
}

void QucsDoc::slotToolBar(QAction* a)
{
	incomplete();
}

// really?!
void QucsDoc::undo()
{ untested();
	QUndoStack* u = undoStack();
	if(u){itested();
		u->undo();
		assert(_app);
		assert(_app->redo);
		_app->redo->setEnabled(true); // yikes.
	}else{ untested();
	}
}

void QucsDoc::redo()
{ untested();
	QUndoStack* u = undoStack();
	if(u){ untested();
		u->redo();
	}else{ untested();
	}
}

QString QucsDoc::fileSuffix (const QString& Name)
{ untested();
  QFileInfo Info (Name);
  return Info.suffix();
}

QString QucsDoc::fileSuffix (void)
{ untested();
  return fileSuffix (_name);
}

QString QucsDoc::fileBase (const QString& Name)
{ untested();
  QFileInfo Info (Name);
  return Info.completeBaseName();
}

QString QucsDoc::fileBase (void)
{ untested();
  return fileBase (_name);
}

#if 0
// cleanup debris
QAction* QucsDoc::selectAction()
{ untested();
	unreachable();
	assert(_app);
  	return _app->select;
}
#endif

// set mode??
void QucsDoc::setActiveAction(MouseAction* a)
{ untested();
	if(eventHandler()){ untested();
		eventHandler()->setActive(a);
	}else{ untested();
	}
}

#if 0
MouseActions const* QucsDoc::mouseActions() const
{ untested();
	auto ma = const_cast<MouseActions*>(mouseActions());
	return ma;
}

MouseActions* QucsDoc::mouseActions()
{ untested();
	return nullptr;
}
#endif

QucsApp* QucsDoc::app()
{ untested();
	QObject* w = dynamic_cast<QWidget*>(this);
	while(w){ untested();
		if(auto a = dynamic_cast<QucsApp*>(w)){ untested();
			return a;
		}else{ untested();
			w = w->parent();
		}
	}
}

// QWidget* QucsDoc::workToolbar()
// { untested();
// 	if(app()){ untested();
// 		return app()->workToolbar();
// 	}else{ untested();
// 		return nullptr;
// 	}
// }
void QucsDoc::printCursorPosition(int x, int y)
{
	if(app()){ untested();
	// BUG. signal status bar?
		app()->printCursorPosition(x, y);
	}else{ untested();
	}
}

void QucsDoc::showEvent(QShowEvent* ev)
{
	if(_eventHandler){
		_eventHandler->setControls(this);
	}else{
	}
}

void QucsDoc::hideEvent(QHideEvent* ev)
{ untested();
	if(app()){ untested();
		// necessary?
		// when hiding a qucsdoc, shouldn't all children become invisible?
		// (does not seem to work)
		app()->clearWorkToolbar();
	}else{ untested();
	}
}

// void QucsDoc::addToolBar(QToolBar* a)
// {
// }

void QucsDoc::addToolBarAction(QAction* a)
{ untested();
	assert(a);
	if(auto w=dynamic_cast<QWidget*>(this)){ untested();
		// assert(a->parent()==w);
		// no. move dto schematicAction
	}else{ untested();
	}

	if(app()){ untested();
		// TODO: add toolbar in one piece.
		app()->addWorkToolbarAction(a);
	}else{ untested();
	}

	if(auto w=dynamic_cast<QWidget*>(this)){ untested();
		// should not have messed with it.
		// assert(a->parent()==w);
	}else{ untested();
	}
}
// similar to former Schematic::performtoggleaction. but take care of actions,
// and deal with undoable commands.
// https://www.walletfox.com/course/qundocommandexample.php?
// https://stackoverflow.com/questions/32973326/qgraphicsscene-how-to-map-item-movements-into-qundocommand
void QucsDoc::possiblyToggleAction(MouseAction* a, QObject* s)
{ untested();
	auto sender = dynamic_cast<QAction*>(s);
	if(!sender){ // assert?
		incomplete();
		return;
	}else{
	}

	trace3("possiblyToggle", a, sender, sender->text());
	QUndoCommand* cmd = nullptr;
	assert(a);
	if(!sender){ untested();
		setActiveAction(nullptr);
		// cmd = a->activate(sender);
		setActiveAction(a);
	}else if(!sender->isCheckable()){ untested();
		cmd = a->activate(sender);
	}else if(sender->isChecked()){itested();
		cmd = a->activate(sender);

		if(cmd){itested();
			sender->blockSignals(true);
			sender->setChecked(false);
			sender->blockSignals(false);
			// possible 'delete' after select.
			// don't do anything else
			a->deactivate();
		}else{itested();
			// sender->setChecked(true); // assert?
			setActiveAction(a);
		}
	}else{itested();
		setActiveAction(nullptr);
	}

	if(cmd){itested();
		executeCommand(cmd);
	}else{itested();
	}
}
/* -------------------------------------------------------------------------------- */
// maybe this only works for SchematicDoc.
// SchematicDoc has input modes coupled to "MouseActions" that deal with user input.
// TextDoc also has modes, but somehow redirects input to another widget.
MouseAction* QucsDoc::activeAction()
{ untested();
	if(eventHandler()){ untested();
		return eventHandler()->activeAction();
	}else{ untested();
		return nullptr;
	}
}
/* -------------------------------------------------------------------------------- */
MouseAction const* QucsDoc::activeAction() const
{ untested();
	auto d = const_cast<QucsDoc*>(this);
	return d->activeAction();
}
/* -------------------------------------------------------------------------------- */
void QucsDoc::addElement(Element*)
{ untested();
	incomplete();
}
/* -------------------------------------------------------------------------------- */
SchematicModel* QucsDoc::model()
{ untested();
	incomplete();
	return nullptr;
}
/* -------------------------------------------------------------------------------- */
void QucsDoc::executeCommand(QUndoCommand* c)
{ untested();
	// TODO: what if there are multiple views to a scene?
	// is mouseActions == scene?

	if(eventHandler()){itested();
		eventHandler()->executeCommand(c); // really?
		// setChanged();
		if(!DocChanged){ untested();
			emit signalFileChanged(true);
		}else{ untested();
			// TODO: unset Changed in undo.
		}
		DocChanged = c;

		showBias = -1;   // schematic changed => bias points may be invalid
	}else{ untested();
	}
}
/* -------------------------------------------------------------------------------- */
CommonData* QucsDoc::qucsData(std::string const& /*key*/)
{ untested();
	return nullptr; // _data[key];
}
/* -------------------------------------------------------------------------------- */
void QucsDoc::slotDCbias()
{ untested();
	incomplete();
}
/* -------------------------------------------------------------------------------- */
void QucsDoc::slotSimulate()
{ untested();
  incomplete();
} // QucsDoc::slotSimulate
/* -------------------------------------------------------------------------------- */
// move around more.
#include <QFileDialog>
#include <QMessageBox>
bool QucsDoc::saveAs()
{ untested();
	int n = -1;
	QString s;
	QString Filter = QWidget::tr("Schematic") + " (*.sch);;" + QWidget::tr("Any File")+" (*)";
	QFileInfo Info;
	while(true) { untested();
		s = docName();
		Info.setFile(s);
		// which is default directory ?
		if(!s.isEmpty()) { untested();
#if 0
		}else if(!ProjName.isEmpty()) { untested();
			s = QucsSettings.QucsWorkDir.path();
		}else if(!lastDirOpenSave.isEmpty()){ untested();
			s = lastDirOpenSave;
#endif
		} else{ untested();
			s = QDir::currentPath();
		}

		// list of known file extensions
		QString ext = "vhdl;vhd;v;va;sch;dpl;m;oct;net;qnet;txt";
		QStringList extlist = ext.split (';');

		auto w = dynamic_cast<QWidget*>(this);
		assert(w);

		s = QFileDialog::getSaveFileName(w, QWidget::tr("Enter a Document Name"),
				QDir(QString::fromStdString(QucsSettings.QucsWorkDir)).absolutePath(),
				Filter);

		if(s.isEmpty()) { untested();
			return false;
		}else{ untested();
		}

		Info.setFile(s);               // try to guess the best extension ...
		ext = Info.suffix();

		if(ext.isEmpty() || !extlist.contains(ext)) { untested();
			incomplete();
			// if no extension was specified or is unknown
		}else{ untested();
		}

		Info.setFile(s);
		if(QFile::exists(s)) { untested();
			n = QMessageBox::warning(w, QWidget::tr("Warning"),
					QWidget::tr("The file '")+Info.fileName()+QWidget::tr("' already exists!\n")+
					QWidget::tr("Saving will overwrite the old one! Continue?"),
					QWidget::tr("No"), QWidget::tr("Yes"), QWidget::tr("Cancel"));
			if(n == 2){ untested();
				return false;    // cancel
			}else if(n == 0) { untested();
				// try again?
				continue;
			}else{ untested();
			}
		}else{ untested();
		}

#if 0 // incomplete
		// search, if document is open
		QucsDoc * d = findDoc (s);
		if(d) { untested();
			QMessageBox::information(this, QWidget::tr("Info"),
					QWidget::tr("Cannot overwrite an open document"));
			return false;
		}
#endif

		break;
	}
	setName(s);

	auto W = dynamic_cast<QWidget*>(this);

	if(auto o=prechecked_cast<QTabWidget*>(_owner)){ untested();
		o->setTabText(o->indexOf(W), s); // misc::properFileName(s));
		// lastDirOpenSave = Info.absolutePath();  // remember last directory and file
	}else{ untested();
	}


	n = save();   // SAVE
	if(n < 0){ untested();
		return false;
	}else{ untested();
	}

//	updatePortNumber(Doc, n);
//
//	// signal?
//	updateRecentFilesList(s);
	return true;
}
/* -------------------------------------------------------------------------------- */
#include "qucs_globals.h"
#include "simulator.h"
/* -------------------------------------------------------------------------------- */
void Simulator::attachDoc(QucsDoc /* const?? */ * d)
{ untested();
	_doc = d;
}
/* -------------------------------------------------------------------------------- */
#if 1
// FIXME.
Simulator* QucsDoc::simulatorInstance(std::string const& which)
{ untested();
	if (!scope()){ untested();
		return nullptr;
	}else{ untested();
	}
	auto i = scope()->find_(which);
	Simulator* sim=nullptr;

	if(i==scope()->end()){ untested();

		Simulator const* proto = QucsSettings.simulator();
		if(which!=""){ untested();
			auto d = dataDispatcher[which];
			proto = dynamic_cast<Simulator*>(d);
		}else{ untested();
		}

		if(proto){ untested();
			sim = proto->clone();
			sim->set_owner(root());
			scope()->push_back(sim);
		}else{ untested();
			incomplete();
		}
	}else if(auto s=dynamic_cast<Simulator*>(*i)){ untested();
		sim = s;
	}else{ untested();
	}

	return sim;
}
#endif
/* -------------------------------------------------------------------------------- */
void setParameter(std::string const&, std::string const&)
{ untested();
	incomplete();
}
/* -------------------------------------------------------------------------------- */
SchematicModel* QucsDoc::scope()
{ untested();
	if(root()){ untested();
		return root()->scope();
	}else{ untested();
		return nullptr;
	}
}
/* -------------------------------------------------------------------------------- */
SchematicModel const* QucsDoc::scope() const
{ untested();
	auto d = const_cast<QucsDoc*>(this);
	return d->scope();
}
/* -------------------------------------------------------------------------------- */
SubcktBase const* QucsDoc::root() const
{ untested();
	auto d = const_cast<QucsDoc*>(this);
	return d->root();
}
/* -------------------------------------------------------------------------------- */
#if 0
// don't know how this works in legacy qucs. the tabs suggest that circuits
// from different "projects" can be openened simultaneously.
// Now a project might have prototypes (e.g. user library) that are not global,
// they are stashed here.
void QucsDoc::installElement(Element const* e)
{ untested();
	trace1("installElement", e->label());
	auto& p = _protos[e->label()];
	if(p){ untested();
		delete p;
	}else{ untested();
	}
	p = e;
}
#endif
/* -------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------- */
