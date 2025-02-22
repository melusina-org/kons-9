(in-package #:kons-9)

;;;; point-cloud ========================================================

(defclass point-cloud (shape)
  ((points :accessor points :initarg :points :initform (make-array 0 :adjustable t :fill-pointer t))
   (point-colors :accessor point-colors :initarg :point-colors :initform nil)))

(defmethod printable-data ((self point-cloud))
  (strcat (call-next-method) (format nil ", ~a points" (length (points self)))))

(defmethod get-bounds ((p-cloud point-cloud))
  (when (= 0 (length (points p-cloud)))
    (warn "POINT-CLOUD ~a does not have any points. Using default bounds values." p-cloud)
    (return-from get-bounds (values (p! -1 -1 -1) (p! 1 1 1))))
  (points-bounds (points p-cloud)))

;;; TODO -- not tested
(defmethod get-global-bounds ((p-cloud point-cloud))
  (when (= 0 (length (points p-cloud)))
    (warn "POINT-CLOUD ~a does not have any points. Using default bounds values." p-cloud)
    (return-from get-global-bounds (values (p! -1 -1 -1) (p! 1 1 1))))
  (let* ((points (points p-cloud))
         (bounds-lo (p:copy (aref points 0)))
         (bounds-hi (p:copy (aref points 0))))
    (if (scene p-cloud)
        (let* ((paths (get-shape-paths (scene p-cloud) p-cloud))
               (path (first paths))
               (matrix (if path (shape-global-matrix (scene p-cloud) path) nil)))
          (if matrix
              (progn
                (do-array (i p points)
                  (let ((xform-p (transform-point p matrix)))
                    (p:min! bounds-lo bounds-lo xform-p)
                    (p:max! bounds-hi bounds-hi xform-p)))
                (values bounds-lo bounds-hi))
              (get-bounds p-cloud)))
        (get-bounds p-cloud))))

(defun make-point-cloud (points &optional (colors nil))
  (make-instance 'point-cloud :points points :point-colors colors))

(defmethod append-point ((p-cloud point-cloud) point &optional (color nil))
  (vector-push-extend point (points p-cloud))
  (when (point-colors p-cloud)
    (vector-push-extend (or color (fg-color *drawing-settings*)) (point-colors p-cloud)))
  p-cloud)

(defmethod freeze-transform ((p-cloud point-cloud))
  (transform-point-array! (points p-cloud) (transform-matrix (transform p-cloud)))
  (reset-transform (transform p-cloud))
  p-cloud)

(defmethod allocate-point-colors ((p-cloud point-cloud)
                                  &optional (color (fg-color *drawing-settings*)))
  (setf (point-colors p-cloud) (make-array (length (points p-cloud))
                                           :adjustable t :fill-pointer t
                                           :initial-element color))
  p-cloud)

(defmethod reset-point-colors ((p-cloud point-cloud))
  (allocate-point-colors p-cloud)
  p-cloud)

(defmethod set-point-colors ((p-cloud point-cloud) color)
  (allocate-point-colors p-cloud color))

(defmethod set-point-colors-by-xyz ((p-cloud point-cloud) color-fn)
  (allocate-point-colors p-cloud)
  (do-array (i p (points p-cloud))
    (setf (aref (point-colors p-cloud) i) (funcall color-fn p)))
  p-cloud)

(defmethod set-point-colors-by-order ((p-cloud point-cloud) color-fn)
  (allocate-point-colors p-cloud)
  (let ((n (length (points p-cloud))))
    (do-array (i p (points p-cloud))
      (declare (ignore p))
      (setf (aref (point-colors p-cloud) i) (funcall color-fn (/ i n)))))
  p-cloud)

;;; point generator functions --------------------------------------------------

(defun make-line-points (p1 p2 num-segments)
  (let* ((num-points (1+ num-segments))
         (points (make-array num-points)))
    (dotimes (i num-points)
      (setf (aref points i) (p:lerp p1 p2 (/ i (coerce num-segments 'single-float)))))
    points))

(defun make-rectangle-points (width height &optional (num-segments 1))
  (let* ((points (make-array (* 4 num-segments)))
         (index -1)
         (x (/ width 2))
         (y (/ height 2))
         (p0 (p!    x     y  0))
         (p1 (p! (- x)    y  0))
         (p2 (p! (- x) (- y) 0))
         (p3 (p!    x  (- y) 0))
         (side-0 (make-line-points p0 p1 num-segments))
         (side-1 (make-line-points p1 p2 num-segments))
         (side-2 (make-line-points p2 p3 num-segments))
         (side-3 (make-line-points p3 p0 num-segments)))
    ;; skip last point of each side
    (dotimes (i num-segments)
      (setf (aref points (incf index)) (aref side-0 i)))
    (dotimes (i num-segments)
      (setf (aref points (incf index)) (aref side-1 i)))
    (dotimes (i num-segments)
      (setf (aref points (incf index)) (aref side-2 i)))
    (dotimes (i num-segments)
      (setf (aref points (incf index)) (aref side-3 i)))
    points))

(defun make-circle-points (diameter num-segments)
  (let ((points (make-array num-segments))
        (radius (/ diameter 2.0))
        (angle-delta (/ 2pi num-segments)))
    (dotimes (i num-segments)
      (let ((angle (* i angle-delta)))
        (setf (aref points i) (p! (* (sin angle) radius) (* (cos angle) radius) 0))))
    (nreverse points)))                 ;return ccw points

(defun make-arc-points (diameter start-angle end-angle num-segments)
  (let* ((points (make-array (1+ num-segments)))
         (radius (/ diameter 2.0))
         (angle-delta (/ (- (radians end-angle) (radians start-angle)) num-segments)))
    (dotimes (i (1+ num-segments))
      (let ((angle (+ (* i angle-delta) (radians start-angle))))
        (setf (aref points i) (p! (* (sin angle) radius) (* (cos angle) radius) 0))))
    points))

(defun make-spiral-points (start-diameter end-diameter axis-length num-loops num-segments)
  (let ((points (make-array (1+ num-segments)))
        (start-radius (/ start-diameter 2.0))
        (end-radius (/ end-diameter 2.0))
        (angle-delta (/ (* 2pi num-loops) num-segments))
        (len-delta (/ axis-length num-segments)))
    (dotimes (i (1+ num-segments))
      (let ((angle (* i angle-delta))
            (len (* i len-delta))
            (radius (lerp (/ i num-segments) start-radius end-radius)))
        (setf (aref points i) (p! (* (sin angle) radius) (* (cos angle) radius) len))))
    (nreverse points)))                 ;return ccw points

(defun make-sine-curve-points (period frequency x-scale y-scale num-segments)
  (let* ((points (make-array (1+ num-segments)))
         (rad-period (radians period))
         (angle-delta (/ rad-period num-segments)))
    (dotimes (i (1+ num-segments))
      (let ((angle (* i angle-delta frequency)))
        (setf (aref points i) (p! (* x-scale (/ angle (* frequency rad-period))) (* y-scale (sin angle)) 0))))
    points))

(defun make-random-points (num bounds-lo bounds-hi)
  (let ((points (make-array num)))
    (dotimes (i num)
      (setf (aref points i) (p-rand2 bounds-lo bounds-hi)))
    points))

(defun make-grid-points (nx ny nz bounds-lo bounds-hi)
  (let ((points (make-array (* nx ny nz)))
        (i -1))
    (dotimes (ix nx)
      (let* ((fx (tween ix 0 (- nx 1.0)))
	     (x (lerp fx (p:x bounds-lo) (p:x bounds-hi))))
        (dotimes (iy ny)
          (let* ((fy (tween iy 0 (- ny 1.0)))
                 (y (lerp fy (p:y bounds-lo) (p:y bounds-hi))))
            (dotimes (iz nz)
              (let* ((fz (tween iz 0 (- nz 1.0)))
                     (z (lerp fz (p:z bounds-lo) (p:z bounds-hi))))
                (setf (aref points (incf i)) (p! x y z))))))))
    points))

(defmethod randomize-points ((p-cloud point-cloud) delta)
  (do-array (i p (points p-cloud))
    (let ((offset (p! (rand1 (p:x delta)) (rand1 (p:y delta)) (rand1 (p:z delta)))))
      (p:+! p p offset)))
  p-cloud)
