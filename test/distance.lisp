(in-package #:kons-9)

(defparameter *tetrahedron*
  (translate-to (make-tetrahedron  2.0) (p! -5 0 0)))

(defparameter *cube*
  (translate-to (make-cube         2.0) (p! -2.5 0 0)))

(defparameter *octahedron*
  (translate-to (make-octahedron   2.0) (p! 0 0  0)))

(defparameter *dodecahedron*
  (translate-to (make-dodecahedron 2.0) (p!  2.5 0 0)))

(defparameter *icosahedron*
  (translate-to (make-icosahedron  2.0) (p! 5 0 0)))

(defparameter *polyhedrons*
  (list
   *tetrahedron*
   *cube*
   *octahedron*
   *dodecahedron*
   *icosahedron*))

(defun polyhedron->coordinate-matrix (polyhedron)
  (let* ((n
	   (length (points polyhedron)))
	 (matrix
	   (make-array (list 3 n) :initial-element 0.0)))
    (loop :for point :across (points polyhedron)
	  :for i = 0 :then (1+ i)
	  :do (setf (aref matrix 0 i) (x point)
		    (aref matrix 1 i) (y point)
		    (aref matrix 2 i) (z point)))
    matrix))

(defun ones (designator)
  (etypecase designator
    (integer
     (make-array designator :initial-element 1.0))
    (array
     (make-array (array-dimensions designator) :initial-element 1.0))
    (list
     (mapcar (lambda (_) (declare (ignore _)) 1.0) designator))))

(defun p-weighted-center (points weights)
  (loop :with weighted-center = +origin+
	:for point :across points
	:for weight :across weights
	:for accumulated-weight = weight :then (+ accumulated-weight weight)
	:do (setf weighted-center
		  (p+ weighted-center
		      (p* point weight)))
	:finally (return (p/ weighted-center accumulated-weight))))

(defun polyhedron-distance (points-1 points-2)
  (let* ((size-1
	   (length points-1))
	 (size-2
	   (length points-2))
	 (weights-1
	   (make-array size-1 :initial-element (/ 1.0 size-1)))
	 (weights-2
	   (make-array size-2 :initial-element (/ 1.0 size-2)))
	 (gradient-1
	   (make-array size-1 :initial-element 0.0))
	 (gradient-2
	   (make-array size-2 :initial-element 0.0))
	 (delta
	   +origin+)
	 (epsilon
	   +origin+))
    (labels ((project-vector-onto-standard-simplex (vector size)
	       (let ((height-over-size
		       (/ (loop :for coordinate :across vector
				:sum coordinate)
			  size)))
		 (loop :for i :from 0 :to (1- size)
		       :do (decf (aref vector i) height-over-size))))
	     (compute-delta ()
	       (setf delta
		     (p- (p-weighted-center points-1 weights-1)
			 (p-weighted-center points-2 weights-2))))
	     (compute-gradient-1 ()
	       (dotimes (i size-1)
		 (setf (aref gradient-1 i)
		       (p-dot (aref points-1 i) delta))))
	     (compute-gradient-2 ()
	       (dotimes (i size-2)
		 (setf (aref gradient-2 i)
		       (- (p-dot (aref points-2 i) delta)))))
	     (compute-epsilon ()
	       (setf epsilon
		     (p- (p-weighted-center points-1 gradient-1)
			 (p-weighted-center points-2 gradient-2))))
	     (minimize-q (delta epsilon)
	       (let ((t-optimum
		       (- (/ (p-dot delta epsilon) 
			     (p-dot epsilon epsilon))))
		     (t-min
		       0.0)
		     (t-max
		       0.0)))))
      (compute-delta)
      (compute-gradient-1)
      (compute-gradient-2)
      (compute-epsilon)
      (format t "~&DELTA ~A~&EPSILON ~A~&GRADIENT-1 ~A~&GRADIENT-2 ~A"
	      delta epsilon gradient-1 gradient-2))))
	       
	       


(p-weighted-center (points *tetrahedron*) (ones (points *tetrahedron*)))

(dolist (polyhedron *polyhedrons*)
  (polyhedron-bake polyhedron))

(with-clear-scene
  (add-shapes *scene* *polyhedrons*))
