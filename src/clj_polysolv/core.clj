(ns clj-polysolv.core
  (:require [clojure.math.numeric-tower :as math]
            [clojure.set :as set]))

; math

(defn round [n precision]
  (let [factor (math/expt 10 precision)]
    (double (/ (Math/round (double (* factor n)))
               factor))))

(defn average [coll]
  (/ (reduce + coll)
     (count coll)))

; polynomial related functions

(defn highest-power [f]
  (if (pos? (count f))
      (apply max (keys f))
      0))

(defn normalize [f]
  (into {}
    (for [i (range (inc (highest-power f)))]
      [i (get f i 0.0)])))

(defn derivative [f]
  (into {}
    (for [[i value] (dissoc f 0)]
      [(dec i) (* i value)])))

; calculate value of polynomial for a given x
(defn poly [f x]
  (reduce +
    (for [[i value] f]
      (* value (math/expt x i)))))

; root solving functions

(defn linear-roots [f]
  (let [f (normalize f)
        a (f 1)
        b (f 0)]
    #{(/ (- b) a)}))

(defn discriminant [a b c]
  (- (math/expt b 2) (* 4 a c)))

(defn quadratic-root [sign a b disc]
  (/ (sign (- b) (math/sqrt disc)) (* 2 a)))

(defn quadratic-roots [f]
  (let [f (normalize f)
        a (f 2)
        b (f 1)
        c (f 0)
        disc (discriminant a b c)]
    (cond
      (pos? disc)  #{(quadratic-root + a b disc)
                     (quadratic-root - a b disc)}
      (zero? disc) #{(quadratic-root + a b disc)}
      (neg? disc)  #{})))

; newton's method guessing functions

(defn good-enough? [prev guess]
  (= (round prev 10) (round guess 10)))

(defn improve-guess [f guess]
  (let [f' (derivative f)
        x  (poly f guess)
        x' (poly f' guess)]
    (if (zero? x')
      :zero-derivative
      (- guess (/ x x')))))

(defn newton-guess [f init]
  (loop [prev init
         guess (improve-guess f init)]
    (cond
      ; guess converged to a division by zero
      ; this is an error condition
      ; backtrack and try a different initial value
      (= :zero-derivative guess)
        (recur f (+ init 1.0))
      (good-enough? prev guess)
        (round guess 4)
      :else
        (recur guess (improve-guess f guess)))))

; helpers

(defn sign-same? [a b]
  (= (pos? a) (pos? b)))

(defn sign-change? [a b]
  (not= (pos? a) (pos? b)))

(defn flatten-start-values [xs]
  (->> xs
       flatten
       (filter identity)
       set))

(defn detect-start-values [f extrema]
  (let [power (highest-power f)
        a (f power)
        first-extremum (first extrema)
        last-extremum (last extrema)]

    (flatten-start-values
      (list
        ; search to the left of first extremum
        ;    positive ^3 AND extremum > 0
        ; OR negative ^3 AND extremum < 0
        ; OR positive ^4 AND extremum < 0
        ; OR negative ^4 AND extremum > 0
        (let [sign-check (if (even? power) sign-change? sign-same?)]
          (when (sign-check a (poly f first-extremum))
            (- first-extremum 1)))

        ; search in between extrema
        ; pairwise iteration
        (map (fn [x1 x2]
                (when (or (zero? (poly f x1))
                          (zero? (poly f x2))
                          (sign-change? (poly f x1) (poly f x2)))
                  (average [x1 x2])))
             (drop-last extrema)
             (rest extrema))

        ; search to the right of last extremum
        ;    positive ^3 AND extremum < 0
        ; OR positive ^4 AND extremum < 0
        ; OR negative ^3 AND extremum > 0
        ; OR negative ^4 AND extremum > 0
        (when (sign-change? a (poly f last-extremum))
          (+ last-extremum 1))))))

(declare solve)

; newton's method (using gauss–lucas theorem)
;
; recursively searches for polynomial roots
;
; y = ax^3 + bx^2 + cx + d
; num roots possible 1-3
;
; y = ax^n + bx^(n-1) + cx^(n-2) ... zx^0
; num roots possible
;        if n even: 0-n
;        if n odd: 1-n
(defn newton-roots [f]
  (let [f (normalize f)
        power (highest-power f)
        f' (derivative f)
        extrema (apply sorted-set (solve f'))]

    ; of odd power with no extrema
    ; only one solution possible
    ; search for it
    (if (and (odd? power) (zero? (count extrema)))
        #{(newton-guess f 1.0)}
        (set (map (partial newton-guess f)
                  (detect-start-values f extrema))))))

(defn solve [f]
  (let [power (highest-power f)]
    (cond
      (< power 2) (linear-roots f)
      (< power 3) (quadratic-roots f)
      :else (newton-roots f))))
