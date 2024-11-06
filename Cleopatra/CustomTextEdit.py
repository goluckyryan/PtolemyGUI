#!/usr/bin/python3

from PyQt6.QtWidgets import QTextEdit
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont, QTextCharFormat, QSyntaxHighlighter

class PythonHighlighter(QSyntaxHighlighter):
  def __init__(self, document):
    super().__init__(document)

    # Define formatting for comments
    self.comment_format = QTextCharFormat()
    self.comment_format.setForeground(Qt.GlobalColor.darkGreen)

  def highlightBlock(self, text):
    # Highlight comments
    if text.startswith("#"):
      self.setFormat(0, len(text), self.comment_format)
    if text.startswith("$"):
      self.setFormat(0, len(text), self.comment_format)
    if text.startswith("0"):
      self.setFormat(0, len(text), self.comment_format)

class CustomTextEdit(QTextEdit):
  def __init__(self, parent=None):
    super().__init__(parent)

    self.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
    font = QFont("Courier New", 10)  # You can adjust the size as needed
    self.setFont(font)

    self.highlighter = PythonHighlighter(self.document())

  def keyPressEvent(self, event):
    # Check if Ctrl+D is pressed
    if event.key() == Qt.Key.Key_D and event.modifiers() == Qt.KeyboardModifier.ControlModifier:
      self.duplicate_line()
      event.accept()  # Prevent the default behavior of Ctrl+D
    else:
        super().keyPressEvent(event)  # Call the base class to handle other keys

  def duplicate_line(self):
    cursor = self.textCursor()

    # Select the current line under the cursor
    cursor.select(cursor.SelectionType.LineUnderCursor)
    line_text = cursor.selectedText()

    # Move the cursor to the end of the line and insert a newline with duplicated text
    cursor.movePosition(cursor.MoveOperation.EndOfLine)
    cursor.insertText("\n" + line_text)
    
    # Update the cursor position after inserting the text
    self.setTextCursor(cursor)
