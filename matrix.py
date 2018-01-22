#
# Implementation of a sparse matrix
#

class Node(object):
  """Objects of type Node represent all non-zero entries of the matrix.
A Node object stores the coordinates of the entry, its value el,
and has a link to the next non-zero entry to the right (in the same row)
and below (in the same column)."""
  def __init__(self, row, col, el, right, down):
    self.row = row
    self.col = col
    self.el = el
    self.right = right
    self.down = down

class Matrix(object):
  def __init__(self, nrows, ncols):
    self.nrows = nrows
    self.ncols = ncols
    self._prow = [None] * nrows
    self._pcol = [None] * ncols

  def _findnode(self, row, col):
    """Returns the node for (row, col) and the previous node in the same row.
Both are None if they do not exist."""
    p = self._prow[row]
    q = None
    while p is not None and p.col < col:
      q = p
      p = p.right
    if p is None or p.col == col:
      return p, q
    return None, q

  def _insertnode(self, row, col, q, el):
    """Insert a new node for entry (row, col) with value el.
q is the previous node on the same row, or None."""
    n = Node(row, col, el, None, None)
    if q is not None:
      n.right = q.right
      q.right = n
    else:
      n.right = self._prow[row]
      self._prow[row] = n
    p = self._pcol[col]
    m = None
    while p is not None and p.row < row:
      m = p
      p = p.down
    if p is None or p.row == row:
      if m is not None:
        m.down = n
      else:
        self._pcol[col] = n
    else:
      n.down = p
      if m is None:
        self._pcol[col] = n
      else:
        m.down = n

  def _removenode(self, p, q):
    "Remove the node p. q is the previous node on the same row, or None."
    if q is not None:
      if p.right is not None:
        q.right = p.right
      else:
        q.right = None
    else:
      if p.right is not None:
        self._prow[p.row] = p.right
      else:
        self._prow[p.row] = None
    if self._pcol[p.col] == p:
      if p.down is not None:
        self._pcol[p.col] = p.down
      else:
        self._pcol[p.col] = None
    else:
      n = self._pcol[p.col]
      m = None
      while n.row != p.row:
        m = n
        n = n.down
      if p.down is not None:
        m.down = p.down
      else:
        m.down = None

  def __getitem__(self, pos):
    "Return matrix entry pos = (row, col)."
    row, col = pos
    p, q = self._findnode(row, col)
    if p is None:
      return 0.0
    return p.el

  def __setitem__(self, pos, el):
    "Set matrix entry pos = (row, col) to value el."
    row, col = pos
    p, q = self._findnode(row, col)
    if p is None:
      if el != 0.0:
        self._insertnode(row, col, q, el)
    else:
      if el == 0.0:
        self._removenode(p, q)
      else:
        p.el = el
    
  def __repr__(self):
    s = ""
    for row in range(min(self.nrows, 10)):
      if row == 0:
        s += "/"
      elif row == self.nrows-1:
        s += "\\"
      else:
        s += "|"
      for col in range(min(self.ncols, 10)):
        s += "%6s " % self[row, col]
      if self.ncols > 10:
        s += "... "
      if row == 0:
        s += "\\\n"
      elif row == self.nrows-1:
        s += "/\n"
      else:
        s += "|\n"
    if self.nrows > 10:
      s += "...\n"
    return s

  def __eq__(self, rhs):
    "Test two matrices for equality."
    if self.nrows != rhs.nrows or self.ncols != rhs.ncols:
      return False
    for row in range(self.nrows):
      p1 = self._prow[row]
      p2 = rhs._prow[row]
      while p1 is not None and p2 is not None:
        if p1.col != p2.col or p1.el != p2.el:
          return False
        p1 = p1.right
        p2 = p2.right
      if p1 is not None or p2 is not None:
        return False
    return True

  def __mul__(self, rhs):
    "Multiply matrix with vector from the right."
    if self.ncols != len(rhs):
      raise ValueError("Dimensions of matrix and vector do not match")
    result = [0.0] * self.nrows
    p = self._prow[0]
    i = 0
    while p is None:
      i += 1
      p = self._prow[i]
    while i != self.nrows-1:
      entry = 0
      q = p
      while q is not None and q.col != self.ncols:
        entry += q.el * rhs[q.col]
        q = q.right
      result[i] += entry
      i += 1
      p = self._prow[i]
    entry = 0
    while p is not None and p.col != self.ncols:
      entry += p.el * rhs[p.col]
      p = p.right
    result[i] += entry
    return result

  def __rmul__(self, lhs):
    "Multiply matrix with vector from the left."
    if self.nrows != len(lhs):
      raise ValueError("Dimensions of matrix and vector do not match")
    result = [0.0] * self.ncols
    p = self._pcol[0]
    i = 0
    while p is None:
      i += 1
      p = self._pcol[i]
    while i != self.ncols-1:
      entry = 0
      q = p
      while q is not None and q.row != self.nrows:
        entry += lhs[q.row] * q.el
        q = q.down
      result[i] += entry
      i += 1
      p = self._pcol[i]
    entry = 0
    while p is not None and p.row != self.nrows:
      entry += lhs[p.row] * p.el
      p = p.down
    result[i] += entry
    return result

  def transposed(self):
    result = Matrix(self.ncols, self.nrows)
    p = self._pcol[0]
    i = 0
    while p is None:
      i += 1
      p = self._pcol[i]
    m = result._prow[i]
    while i != self.ncols - 1:
      result._insertnode(i, p.row, m, p.el)
      q = p.down
      m = result._prow[i]
      while q is not None and q.row != self.nrows:
        result._insertnode(i, q.row, m, q.el)
        m = m.right
        q = q.down
      i += 1
      while self._pcol[i] is None and i != self.ncols - 1:
        i += 1
      p = self._pcol[i]
      m = result._prow[i]
    while p is not None and p.row != self.nrows:
      result._insertnode(i, p.row, m, p.el)
      if m is None:
        m = result._prow[i]
      else:
        m = m.right
      p = p.down    
    return result

  def __add__(self, rhs):
    if self.nrows != rhs.nrows or self.ncols != rhs.ncols:
      raise ValueError("Dimensions of matrices do not match")
    result = Matrix(self.nrows, self.ncols)
    i = 0
    p = self._prow[i]
    while p is None:
      i += 1
      p = self._prow[i]
    m = result._prow[i]
    while i != self.nrows - 1:
      result._insertnode(i, p.col, m, p.el)
      m = result._prow[i]
      q = p.right
      m = result._prow[i]
      while q is not None and q.col != self.ncols:
        result._insertnode(i, q.col, m, q.el)
        m = m.right
        q = q.right
      i += 1
      while self._prow[i] is None and i != self.nrows-1:
        i += 1
      p = self._prow[i]
      m = result._prow[i]
    while p is not None and p.col != self.ncols:
      result._insertnode(i, p.col, m, p.el)
      if m is None:
        m = result._prow[i]
      else:
        m = m.right
      p = p.right
    i = 0
    p = rhs._prow[i]
    while p is None:
      i += 1
      p = rhs._prow[i]
    m = result._prow[i]
    while i != self.nrows - 1:
      q = p
      while q is not None and q.col != self.ncols:
        if m is None:
          result._insertnode(i, q.col, m, q.el)
          m = result._prow[i]
        elif m.col > q.col:
          result._insertnode(i, q.col, None, q.el)
          m = result._prow[i]
        elif m.col == q.col:
          if m.el + q.el == 0:
            a, b = result._findnode(i, m.col)
            result._removenode(a, b)
          else:
            m.el += q.el
        else:
          n = m
          while n is not None and q.col > n.col:
            n = n.right
            if n is None or n.col > q.col:
              result._insertnode(i, q.col, m, q.el)
              m = m.right
              n = m
            elif n.col == q.col:
              if n.el + q.el == 0:
                a, b = result._findnode(i, n.col)
                result._removenode(a, b)
              else:
                n.el += q.el
            if n is not None:
              m = n
        q = q.right        
      i += 1
      while rhs._prow[i] is None and i != self.nrows-1:
        i += 1
      p = rhs._prow[i]
      m = result._prow[i]
    while p is not None and p.col != self.ncols:
      if m is None:
        result._insertnode(i, p.col, m, p.el)
        m = result._prow[i]
      elif m.col > p.col:
        result._insertnode(i, p.col, None, p.el)
        m = result._prow[i]
      elif m.col == p.col:
        if m.el + p.el == 0:
          a, b = result._findnode(i, m.col)
          result._removenode(a, b)
        else:
          m.el += p.el
      else:
        n = m
        while n is not None and p.col > n.col:
          n = n.right
          if n is None or n.col > p.col:
            result._insertnode(i, p.col, m, p.el)
            m = m.right
            n = m
          elif n.col == p.col:
            if n.el + p.el == 0:
              a, b = result._findnode(i, n.col)
              result._removenode(a, b)
            else:
              n.el += p.el
          if n is not None:
            m = n
      p = p.right
    return result
# --------------------------------------------------------------------

def identity(n):
  "Create an nxn identity matrix."
  M = Matrix(n, n)
  for i in range(n):
    M[i,i] = 1.0
  return M

# --------------------------------------------------------------------

if __name__ == "__main__":
  m = identity(4)
  print(m)
  m[1,1] = 7
  print(m)
  m[2,1] = 13
  print(m)
  m[0,3] = -2
  print(m)
  m[3,3] = 0
  print(m)
  m[0,0] = 0
  print(m)
  m2 = Matrix(4, 4)
  m2[0,3] = -2
  m2[1,1] = 7
  m2[2,1] = 13
  print(m2)
  print(m == m2)
  m2[2,2] = 1
  print(m == m2)
  print(m * [ 1, 2, 3, 4 ] )
  print([1, 2, 3, 4] * m)
  t = m.transposed()
  print(t)
  print([1, 2, 3, 4] * t)
  print(t * [ 1, 2, 3, 4 ] )
  m3 = m + t
  print(m3)
  print(m3 == m3.transposed())
  
# --------------------------------------------------------------------
