

class Matrix
  def initialize(*args)
    if args[1].nil?
      @data = args[0]
      @height = @data.length
      @width = @data[0].length
    else
      @width, @height = args
    end
  end

  attr_reader :width, :height, :data

  # convinient inspection
  def to_s
    @data
      .inspect
      .split('],')
      .join("],\n")
  end

  # iterates over each row of matrix
  def each_row
    @data.each_index do |index|
      yield @data[index], index
    end
  end

  # iterates over each column of matrix
  # warning: the complexity is O(n^2)
  def each_col
    0.upto(@width - 1) do |col_number|
      col = @data.collect { |row| row[col_number] }
      yield col col_number
    end
  end

  # addition
  def +(right_operand)
    if @width != right_operand.width ||
      @height != right_operand.height
      raise 'Matrix must have same size'
    end

    additive_operation(right_operand) do |left, right|
      left + right
    end
  end

  # subtraction
  def -(right_operand)
    if @width != right_operand.width ||
      @height != right_operand.height
      raise 'Matrix must have same size'
    end

    additive_operation(right_operand) do |left, right|
      left - right
    end
  end

  # multiplication
  def *(right_operand)
    if @width != right_operand.height
      raise 'Wrong operand'
    end

    new_data = []
    # iterates over the rows of the left matrixt
    0.upto(@height - 1) do |row_index|
      new_data[row_index] = []
      # and all the columns of the right matrix
      0.upto(right_operand.width - 1) do |col_index|
        cell_value = 0
        # and calc scalar product (row * column)
        0.upto(@width - 1) do |index|
          left_value = @data[row_index][index]
          right_value = right_operand.data[index][col_index]
          cell_value += left_value * right_value
        end

        new_data[row_index].push cell_value
      end
    end

    Matrix.new new_data
  end

  def /(right_operand)
    new_data = @data.dup

    if right_operand.is_a?(Numeric)
      0.upto(@width - 1) do |row|
        0.upto(@height - 1) do |col|
          new_data[row][col] /= right_operand
        end
      end
    end

    Matrix.new new_data
  end

  # gaussian elimination
  def triangular_form(other=nil)
    new_data = @data.dup

    # we can really transform only a maximum-size square
    # submatrix, which width = height = min(@width, @height)
    # so...
    tr_size = [@width, @height].min

    0.upto(tr_size - 1) do |col_index|
      triangular_form_partial new_data, col_index, other, false
      triangular_form_partial new_data, col_index, other, true
    end

    Matrix.new new_data
  end

  # calc determinant through the gaussian elimination
  def determinant
    if @width != @height
      raise 'determinant is defined only for the square matricies'
    end

    triangular = triangular_form

    Matrix.determinant_tr(triangular)
  end

  # calc inverse of this matrix
  # Gauss–Jordan elimination
  def inverse
    result = Matrix.identity(@width)

    triangular = triangular_form(result)

    0.upto(@width - 1) do |index|
      if triangular.data[index][index] != 1
        factor = 1.0 / triangular.data[index][index]
        result.data[index].map! { |e| e * factor }
      end
    end

    result
  end

  def self.identity(dim)
    new_data = Array.new(dim) { Array.new(dim, 0) }
    0.upto(dim - 1) do |index|
      new_data[index][index] = 1
    end

    Matrix.new new_data
  end

  private

  def triangular_form_partial(new_data, col_index, other=nil, top=false)
    # first we need to find row with non-zero
    # first element, for the pivot
    if top
      pivot_index = nonzero_row_above(new_data, col_index)
    else
      pivot_index = nonzero_row_below(new_data, col_index)
    end
    # It may not have non-zero element
    # Then all the elements below the main diagonal are zeros,
    # so we're done here
    if top
      cond = pivot_index >= 0
    else
      cond = pivot_index < @height
    end

    if cond
      # otherwise, pivot_index is our guy.
      # after that, the main diagonal element is non-zero
      swap_rows(new_data, pivot_index, col_index)
      swap_rows(other.data, pivot_index, col_index) if other
      pivot_row = new_data[col_index]
      other_pivot_row = other.data[col_index] if other
      pivot_el = pivot_row[col_index].to_f

      # We're moving along the diagonal, so we should start counting
      # index from the number of the current column
      # (but iteration goes over the rows)
      if top
        range = 0.upto(col_index - 1)
      else
        range = (col_index + 1).upto(@height - 1)
      end

      range.each do |row_index|
        next if row_index == col_index

        current_row = new_data[row_index]
        other_current_row = other.data[row_index] if other

        # already zero
        next if current_row[col_index] == 0

        factor = current_row[col_index] / pivot_el
        subtr = pivot_row.collect { |e| e * factor }
        other_subtr = other_pivot_row.collect { |e| e * factor } if other
        # subtract new row from the bottom one
        current_row.map!.with_index do |e, i|
          e - subtr[i]
        end

        if other
          other_current_row.map!.with_index do |e, i|
            e - other_subtr[i]
          end
        end
      end
    end
  end

  # calc determinant of triangular matrix
  def self.determinant_tr(triangular)
    result = 1
    # get mult of all the main diagonal elements
    0.upto(triangular.width - 1) do |index|
      result *= triangular.data[index][index]
    end

    result
  end

  # abstract additive operation. May represent addition or
  # subtraction depending on the yield
  def additive_operation(right_operand)
    new_data = []

    # iterating all over the elements and DO something with them
    each_row do |row, row_index|
      new_row = []
      row.each_index do |col_index|
        left_value = row[col_index]
        right_value = right_operand.data[row_index][col_index]
        new_row.push(yield left_value, right_value)
      end

      new_data.push new_row
    end

    Matrix.new new_data
  end

  def swap_rows(data, i, j)
    data[i], data[j] = data[j], data[i]
  end

  # returns index of first non-zero element in exact row
  # index >= data.length means all the elements are zeros
  def nonzero_row_below(data, col_index)
    offset = 0
    while (col_index + offset < data.length &&
      data[col_index + offset][col_index] == 0)
      offset += 1
    end
    col_index + offset
  end

  # returns index of first non-zero element in exact row
  # index < 0 means all the elements are zeros
  def nonzero_row_above(data, col_index)
    offset = 0
    while (col_index - offset >= 0 &&
      data[col_index - offset][col_index] == 0)
      offset -= 1
    end
    col_index - offset
  end
end

m1 = Matrix.new([[1, 2, 3], [2, 3, 4], [4, 55, 62]])

puts m1.inverse
