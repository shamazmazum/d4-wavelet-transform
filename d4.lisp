(in-package :d4-wavelet-transform)

(defparameter *coeff* #(0.48296294 0.8365163 0.22414388 -0.12940952)) ; c0 c1 c2 c3
(defparameter *coeff2* #(0.12940952 0.22414388 -0.8365163 0.48296294)); -c3 c2 -c1 c0
(defparameter *beginning* #(0.0 1.0 0.0))

;; Visual representation of the scaling function and the mother wavelet
(defun larger-scale (array)
  (let ((len (length array)))
    (if (not (zerop (rem len 3)))
        (error "len/3 is not an integer"))
    (let ((new-array (make-array (* len 2) :initial-element 0)))
      (loop
         for idx from 0 by (floor len 3)
         for j below 4 do
           (loop for i below len do
                (incf (aref new-array (+ i idx))
                      (* (sqrt 2.0) (aref *coeff* j) (aref array i)))))
      new-array)))

(defun n-scale (n &optional (array *beginning*))
  "This function computes scaling function phi(x) given phi (2*x).
Width of support for phi(x) is doubled and you get an array with twice as
many elements as in the original array."
  (if (zerop n) array (n-scale (1- n) (larger-scale array))))

(defun n-wavelet (n)
  "This function computates mother wavelet psi(x) given phi(2*x)"
  (let* ((scaling (n-scale (1- n)))
         (len (length scaling))
         (wavelet (make-array (* 2 len) :initial-element 0)))
    (loop
       for idx from 0 by (floor len 3)
       for i below 4 do
         (loop for j below len do
              (incf (aref wavelet (+ j idx))
                    (* (sqrt 2.0) (aref *coeff2* i) (aref scaling j)))))
    wavelet))

(defun array2file (filename array)
  "Write an array to a file in the format understood by gnuplot."
  (with-open-file (stream
                   filename
                   :direction :output
                   :if-exists :supersede
                   :if-does-not-exist :create)
    (loop
       for i from 0 by 1
       for x across array do
         (format stream "~f ~f~%" i x))))

;; Actual forward and inverse transforms

(defun prepare-array (array)
  ;; Add two elements to the end of the array, as if input data is periodic
  (let* ((len (length array))
         (new-array (make-array (+ len 2))))
    (loop for i below len do
         (setf (aref new-array i) (aref array i)))
    (setf (aref new-array len)
          (aref new-array 0)
          (aref new-array (1+ len))
          (aref new-array 1))
    new-array))

(defun transform-pass (array)
  ;; This function calculates scalar products <f, phi(x-n)> and <f, psi(x-n)>
  ;; given <f, phi(2*x - n)>
  (assert (zerop (rem (length array) 2)))
  (let* ((new-array (prepare-array array))
         (len (length array))
         (phi (make-array (ash len -1) :initial-element 0))
         (psi (make-array (ash len -1) :initial-element 0)))
    (loop
       for i below len by 2
       for idx from 0 by 1 do
         (loop for j below 4 do
              (incf (aref phi idx)
                    (* (aref *coeff* j) (aref new-array (+ i j))))
              (incf (aref psi idx)
                    (* (aref *coeff2* j) (aref new-array (+ i j))))))
    (values phi psi)))

(defun prepare-array-inverse (array)
  (let* ((len (length array))
         (new-array (make-array (1+ len))))
    (loop for i below len do
         (setf (aref new-array (1+ i))
               (aref array i)))
    (setf (aref new-array 0) (aref array (1- len)))
    new-array))

(defun inverse-transform-pass (phi psi)
  ;; This function does just the opposite of transform-pass, given <f, phi(x-n)>
  ;; and <f, psi(x-n)> it returns <f, phi(2*x-n)>
  (assert (= (length phi) (length psi)))
  (let* ((len (length phi))
         (new-phi (prepare-array-inverse phi))
         (new-psi (prepare-array-inverse psi))
         (result (make-array (ash len 1))))
    (loop
       for i below (ash len 1) by 2
       for j from 0 by 1 do
         (setf (aref result i)
               (+ (* (aref *coeff* 2) (aref new-phi j))
                  (* (aref *coeff2* 2) (aref new-psi j))
                  (* (aref *coeff* 0) (aref new-phi (1+ j)))
                  (* (aref *coeff2* 0) (aref new-psi (1+ j))))

               (aref result (1+ i))
               (+ (* (aref *coeff* 3) (aref new-phi j))
                  (* (aref *coeff2* 3) (aref new-psi j))
                  (* (aref *coeff* 1) (aref new-phi (1+ j)))
                  (* (aref *coeff2* 1) (aref new-psi (1+ j))))))
    result))

(defun transform (array &optional result)
  "Perform a d4 wavelet transform"
  (if (= (length array) 1) (push array result)
      (multiple-value-bind (phi psi)
          (transform-pass array)
        (transform phi (push psi result)))))

(defun inverse-transform (coeff)
  "Restore original signal from d4 wavelet transform"
  (if (= (length coeff) 1) (first coeff)
      (inverse-transform
       (cons (inverse-transform-pass (first coeff)
                                     (second coeff))
             (subseq coeff 2)))))

(defun generate-sin (n)
  "Return a sine function, sampled for every x = k * n/(2pi), k = 0...n"
  (let ((delta (/ (* 2 pi) n))
        (array (make-array n)))
    (loop
       for i below n
       for x from 0 by delta do
         (setf (aref array i) (sin x)))
    array))

(defun lose-precision (coeff threshold)
  "Zero coefficients in wavelet transform if their absolute magnitude is less
than threshold"
  (loop for c in coeff do
       (loop for i below (length c) do
            (symbol-macrolet ((x (aref c i)))
              (setf x
                    (if (< (abs x) threshold) 0 x)))))
  coeff)
