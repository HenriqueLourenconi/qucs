/*
 * config-base.h - useful piece of code to include
 *
 * Copyright (C) 2012 Linus Torvarld
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this package; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.  
 *
 * $Id$
 *
 */

#ifndef QUCS_CONFIG_BASE_H
#define QUCS_CONFIG_BASE_H

/* test wether a macro is defined in C
   use like if(is_set(WORDS_BIGENDIAN)) {} else {}
   avoid #ifdef #else #endif
*/
#define is_set(macro) is_set_(macro)
/* internal implementation */
#define macrotest_1 ,
#define is_set_(value) is_set__(macrotest_##value)
#define is_set__(comma) is_set___(comma 1, 0)
#define is_set___(_, v, ...) v


#endif
