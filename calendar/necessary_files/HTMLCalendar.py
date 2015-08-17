#!/usr/local/bin/python

"""HTMLCalendar - Functions and classes for generating one-month and twelve-month calendars in HTML format.

Requires HTMLTemplate.
"""

# HTMLCalendar - Functions and classes for generating one-month and twelve-month calendars in HTML format.
#
# Copyright (C) 2004 HAS <hamish.sanderson@virgin.net>
#
# This library is free software; you can redistribute it and/or modify it under 
# the terms of the GNU Lesser General Public License as published by the Free 
# Software Foundation; either version 2.1 of the License, or (at your option) 
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more 
# details.
#
# You should have received a copy of the GNU Lesser General Public License 
# along with this library; if not, write to the Free Software Foundation, Inc., 
# 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


import calendar, time

from HTMLTemplate import Template


#################################################
# Month Calendar
#################################################

class MonthCal:
	"""One-month HTML calendar renderer."""
	
	_months = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
	_weekdaysFromSun = ["Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"]
	_weekdaysFromMon = _weekdaysFromSun[1:] + [_weekdaysFromSun[0]]
	
	_html = '''<!-- begin one-month calendar -->
<table class="calendar" summary="Monthly calendar.">
<caption node="con:caption"><b>January</b></caption>
<thead>
<tr>
<th node="rep:labels" scope="col" abbr="Sunday">S</th>
<th node="del:" scope="col" abbr="Monday">M</th><th node="del:" scope="col" abbr="Tuesday">T</th><th node="del:" scope="col" abbr="Wednesday">W</th><th node="del:" scope="col" abbr="Thursday">T</th><th node="del:" scope="col" abbr="Friday">F</th><th node="del:" scope="col" abbr="Saturday">S</th>
</tr>
</thead><tbody>
<tr node="rep:week">
<td node="del:">&nbsp;</td><td node="del:">&nbsp;</td>
<td node="rep:day" class="wkday"><a node="con:link" href="archives/2002-01-01.html">1</a></td>
<td class="wkday" node="del:">2</td><td class="wkday" node="del:">3</td><td class="wkday" node="del:">4</td><td class="wkend" node="del:">5</td>
</tr><tr node="del:">
<td class="wkend">6</td>
<td class="wkday">7</td>
<td class="wkday">8</td>
<td class="wkday">9</td>
<td class="wkday">10</td>
<td class="wkday">11</td>
<td class="wkend">12</td>
</tr><tr node="del:">
<td class="wkend">13</td>
<td class="wkday">14</td>
<td class="wkday">15</td>
<td class="wkday">16</td>
<td class="wkday">17</td>
<td class="wkday">18</td>
<td class="wkend">19</td>
</tr><tr node="del:">
<td class="wkend">20</td>
<td class="wkday">21</td>
<td class="wkday">22</td>
<td class="wkday">23</td>
<td class="wkday">24</td>
<td class="wkday">25</td>
<td class="wkend">26</td>
</tr><tr node="del:">
<td class="wkend">27</td>
<td class="wkday">28</td>
<td class="wkday">29</td>
<td class="wkday">30</td>
<td class="wkday">31</td>
<td>&nbsp;</td>
<td>&nbsp;</td>
</tr><tr node="del:">
<td>&nbsp;</td>
<td>&nbsp;</td>
<td>&nbsp;</td>
<td>&nbsp;</td>
<td>&nbsp;</td>
<td>&nbsp;</td>
<td>&nbsp;</td>
</tr>
</tbody>
</table>
<!-- end one-month calendar -->'''

	def __init__(self, isSundayFirst=True):
		"""
			isSundayFirst : boolean -- if True, week begins on Sunday (the default); otherwise Monday
		"""
		self._template = Template(self._renderTemplate, self._html)
		self._isSundayFirst = isSundayFirst
		if isSundayFirst:
			self._firstWeekday = 6
			columnLabels = self._weekdaysFromSun
		else:
			self._firstWeekday = 0
			columnLabels = self._weekdaysFromMon
		self._template.labels.repeat(self._renderLabels, columnLabels) # pre-render column labels for efficiency
		
	def _WeekendTrackerGen(isSundayFirst=True):
		while 1:
			days = [bool(isSundayFirst), False, False, False, False, not isSundayFirst, True]
			while days:
				yield days.pop(0)
	_WeekendTrackerGen = staticmethod(_WeekendTrackerGen)
	
	def _renderLabels(self, node, weekday):
		node.atts['abbr'] = weekday
		node.content = weekday[0]
	
	def _renderTemplate(self, node, year, month, ifo, links, dayrenderer):
		fd = calendar.firstweekday()
		calendar.setfirstweekday(self._firstWeekday)
		weeks = calendar.monthcalendar(year, month) # a list of seven-item lists
		calendar.setfirstweekday(fd)
		if len(weeks) == 5:
			weeks.append([0, 0, 0, 0, 0, 0, 0])
		node.caption.content = "%s %s %s"%(ifo,self._months[month - 1],year) # set table caption
		node.week.repeat(self._renderWeek, weeks, links, self._WeekendTrackerGen(self._isSundayFirst), year, month, dayrenderer) # render weekly rows
	
	def _renderWeek(self, node, days, links, weekendTracker, year, month, dayrenderer):
		node.day.repeat(self._renderDay, days, links, weekendTracker, year, month, dayrenderer)
	
	def _renderDay(self, node, day, links, weekendTracker, year, month, dayrenderer):
		isWeekend = weekendTracker.next()
		if day == 0:
			node.link.raw = '&nbsp;'
			del node.atts['class']
			node.link.omittags()
		else:
			if dayrenderer:
				node.raw = dayrenderer(year, month, day)
			else:
				node.link.content = str(day)
				if isWeekend: 
					node.atts['class'] = 'wkend'
				if links.has_key(day):
					node.link.atts['href'] = links[day]
				else: 
					node.link.omittags()
	
	def render(self, year, month, ifo, links={}, dayrenderer=None):
		# note: the dayrenderer parameter is kinda kludged in; it'd be better if links arg was eliminated and some default day renderers were provided instead, but it's kept for backwards-compatibility's sake
		"""Render a one-month calendar.
			year : integer -- e.g. 2004
			month : integer -- 1-12
			links : dict -- a dict of form {DAY [integer]: URI [string],...} used to add links to selected days; ignored if dayrenderer argument is given; optional
			dayrenderer : function -- a function/closure/callable class instance that takes year, month and day integers as arguments and returns HTML markup for each 'day' table cell; optional
			Result : string -- HTML table
		"""
		return self._template.render(year, month, ifo, links, dayrenderer)



#################################################
# Year Calendar
#################################################

class YearCal:
	"""Twelve-month HTML calendar renderer."""

	_html = '''<!-- begin twelve-month table -->
<table class="yeartable" summary="Year calendar from January to December." width="100%">
<tr node="rep:row">
<td node="rep:col" align="center">MONTH</td>
<td node="rep:col" align="center">MONTH</td>
<td node="rep:col" align="center">MONTH</td>
<td node="rep:col" align="center">MONTH</td>
</tr>
<tr node="rep:row">
<td node="rep:col" align="center">MONTH</td>
<td node="rep:col" align="center">MONTH</td>
<td node="rep:col" align="center">MONTH</td>
<td node="rep:col" align="center">MONTH</td>
</tr>
<tr node="rep:row">
<td node="rep:col" align="center">MONTH</td>
<td node="rep:col" align="center">MONTH</td>
<td node="rep:col" align="center">MONTH</td>
<td node="rep:col" align="center">MONTH</td>
</tr>
</table>
<!-- end twelve-month table -->'''

	def __init__(self, monthCalendar):
		"""
			monthCalendar : MonthCal -- an instance of MonthCal (or any other object that returns a
					one-month HTML calendar table when its render(year, month) method is called)
		"""
		self._template = Template(self._renderTemplate, self._html)
		self._monthCalendar = monthCalendar
	
	def _renderTemplate(self, node, year, dayrenderer):
		node.row.repeat(self._renderRow, [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]], int(year), dayrenderer) # render twelve-month table (4x3)

	def _renderRow(self, node, monthsRow, year, dayrenderer):
		node.col.repeat(self._renderCol, monthsRow, year, dayrenderer)
	
	def _renderCol(self, node, month, year, dayrenderer): # a single month cell
		node.raw = self._monthCalendar.render(year, month, {}, dayrenderer)
	
	def render(self, year=None, dayrenderer=None):
		"""Render a twelve-month calendar as a 4x3 table
			year : integer -- e.g. 2004; default is current year
			dayrenderer : function -- optional dayrenderer function to pass to each month template
			Result : string -- HTML table
		"""
		if year is None:
			year = time.gmtime()[0]
		return self._template.render(year, dayrenderer)


#################################################
# Page renderer
#################################################

_page = '''<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en"><head>

<title node="con:title">%i calendar</title>

<style type="text/css" media="all"><!--
body {padding:0; margin:0;}
h1, .footer {text-align:center; color:#bdb; background-color:#252; padding:4px;}
.calendar {font-weight:bold; color:black; background-color:#484; padding:0px; margin-bottom:12px;}
.calendar caption {font-size:small; font-weight:bold; text-align:center; color:#bdb; background-color:#252; padding:4px;}
.calendar th {text-align:center; color:white; background-color:#252; padding:2px;}
.calendar td {text-align:right; padding:2px;}
.calendar .wkday {color:black; background-color:#bdb;}
.calendar .wkend {color:black; background-color:#9c9;}
--></style>
</head><body>

<h1 node="con:heading">%i</h1>

%s

<div class="footer"><small>Made with HTMLTemplate</small></div>

</body></html>'''

def renderYearCalendar(year=None, isSundayFirst=True):
	"""Render a twelve-month calendar as a complete HTML page.
		year : integer -- e.g. 2004; default is current year
		isSundayFirst : boolean -- if True, week begins on Sunday (the default); otherwise Monday
		Result : string -- HTML document
	"""
	if year is None:
		year = time.gmtime()[0]
	return _page % (year, year, YearCal(MonthCal(isSundayFirst)).render(year))


#################################################
# Test
#################################################

if __name__ == '__main__':
	print renderYearCalendar()

