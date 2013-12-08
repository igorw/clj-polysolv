(ns clj-polysolv.core-test
  (:require [clojure.test :refer :all]
            [clj-polysolv.core :refer :all]))

(deftest derivative-test
  (testing "f = 15x^3 + 10x^2 + 5x
            df = 3*15^2 + 2*10x + 5"
    (let [f {3 15.0,
             2 10.0,
             1 5.0}
          df {2 (* 3 15.0),
              1 (* 2 10.0),
              0 (* 1 5.0)}]
      (is (= df (derivative f))))))

(deftest mock-exam-2a-test
  (testing "-2/3x + 5
            result = 7.5"
    (let [f {1 (/ -2.0 3.0),
             0 5.0}
          result #{7.5}]
      (is (= result (linear-roots f))))))

(deftest mock-exam-2b-test
  (testing "1/2x^2 - 2x - 6
            results = -2, 6"
    (let [f {2 0.5,
             1 -2.0,
             0 -6.0}
          result #{-2.0 6.0}]
      (is (= result (quadratic-roots f))))))

(deftest exam-1a-test
  (testing "-2x^2 + 4x + 6
            results = -1, 3"
    (let [f {2 -2.0,
             1 4.0,
             0 6.0}
          result #{-1.0 3.0}]
      (is (= result (quadratic-roots f))))))

(deftest mock-exam-2e-test
  (testing "x^3 - 3x - 2
            results = -1, 2"
    (let [f {3 1.0,
             1 -3.0,
             0 -2.0}
          result #{-1.0 2.0}]
      (is (= result (newton-roots f))))))

(deftest x-cubed-test
  (testing "x^3
            result = 0"
    (let [f {3 1.0}
          result #{0.0}]
      (is (= result (newton-roots f))))))

(deftest derivative-without-root-test
  (testing "-x^3 - 3x
            result = 0"
    (let [f {3 -1.0,
             1 -3.0}
          result #{0.0}]
      (is (= result (newton-roots f))))))

(deftest nex-x-cubed-test
  (testing "-x^3
            result = 0"
    (let [f {3 -1.0}
          result #{0.0}]
      (is (= result (newton-roots f))))))

(deftest mock-exam-2d-test
  (testing "-1/3x^4 + 4/3x^3
            results = 0, 4"
    (let [f {4 (/ -1 3),
             3 (/ 4 3)}
          result #{0.0 4.0}]
      (is (= result (newton-roots f))))))

(deftest mock-exam-7-test
  (testing "2/3x^3 - 1/2x^2 - 36x + 6
            results = -7.0690, 0.1664, 7.6527"
    (let [f {3 (/ 2 3),
             2 (/ -1 2),
             1 -36,
             0 6}
          result #{-7.0690 0.1664 7.6527}]
      (is (= result (newton-roots f))))))

(deftest mock-exam-2f-test
  (testing "-2.0x^3 + 6x
            results = -1.7320, 0, 1.7320"
    (let [f {3 -2.0,
             1 6.0}
          result #{-1.7321 0.0 1.7321}]
      (is (= result (newton-roots f))))))

(deftest tobias-test
  (testing "1.5x^4 - x^3 - x^2 + 0.1
            results = 0.2930, 1.1880"
    (let [f {4 1.5,
             3 -1.0,
             2 -1.0,
             0 0.1}
          result #{0.2931 1.1881}]
      (is (= result (newton-roots f))))))

(deftest exam-1b-test
  (testing "0.5x^4 - x^3
            results = 0, 2"
    (let [f {4 0.5,
             3 -1.0}
          result #{0.0 2.0}]
      (is (= result (newton-roots f))))))

(deftest exam-1c-test
  (testing "-2x^3 + 3x^2
            results = 0, 1.5"
    (let [f {3 -2.0,
             2 3.0}
          result #{0.0 1.5}]
      (is (= result (newton-roots f))))))

(deftest exam-6-test
  (testing "1/3x^3 + x^2 - 24x + 14
            results = -10.3492, 0.6014, 6.7477"
    (let [f {3 (/ 1 3),
             2 1.0,
             1 -24.0,
             0 14.0}
          result #{-10.3492 0.6014 6.7478}]
      (is (= result (newton-roots f))))))

(deftest high-power-7-test
  (testing "-7/6x^7 + 1/3x^4 + 4/3x^3
            results = -0.96499, 0, 1.09862"
    (let [f {7 (/ -7 6),
             4 (/ 1 3),
             3 (/ 4 3)}
          result #{-0.965 0.0 1.0986}]
      (is (= result (newton-roots f))))))

(deftest derivation-has-only-one-root-test
  (testing "-245x^4 + 8x + 8
            results = -0.377575, 0.467902"
    (let [f {4 -245.0,
             1 8.0,
             0 8.0}
          result #{-0.3776 0.4679}]
      (is (= result (newton-roots f))))))

(deftest no-roots-test
  (testing "53x^4 + 9x^3 + 38x^2 + 4x + 20
            results = none"
    (let [f {4 53.0,
             3 9.0,
             2 38.0,
             1 4.0,
             0 20.0}
          result #{}]
      (is (= result (newton-roots f))))))

(deftest high-power-12-test
  (testing "-6x^12 + 4x^9 - 9x^6 + 2
            results = -0.745162, 0.78614"
    (let [f {12 -6.0,
             9 4.0,
             6 -9.0,
             0 2.0}
          result #{-0.7452 0.7861}]
      (is (= result (newton-roots f))))))

(deftest x-power-4-test
  (testing "x^4
            result = 0"
    (let [f {4 1.0}
          result #{0.0}]
      (is (= result (newton-roots f))))))

(deftest x-power-5-test
  (testing "x^5
            result = 0"
    (let [f {5 1.0}
          result #{0.0}]
      (is (= result (newton-roots f))))))

(deftest x-power-6-test
  (testing "x^6
            result = 0"
    (let [f {6 1.0}
          result #{0.0}]
      (is (= result (newton-roots f))))))

(deftest x-power-7-test
  (testing "x^7
            result = 0"
    (let [f {7 1.0}
          result #{0.0}]
      (is (= result (newton-roots f))))))

(deftest random-0-test
  (testing "x^4 - 81
            results = -3, 3"
    (let [f {4 1.0,
             0 -81.0}
          result #{-3.0 3.0}]
      (is (= result (newton-roots f))))))

(deftest random-1-test
  (testing "x^4 - 6x^2 + 8
            results = +- sqrt(3 +- 1)"
    (let [f {4 1.0,
             2 -6.0,
             0 8.0}
          result #{(round (Math/sqrt 4) 4)
                   (round (Math/sqrt 2) 4)
                   (round (- (Math/sqrt 4)) 4)
                   (round (- (Math/sqrt 2)) 4)}]
      (is (= result (newton-roots f))))))

(deftest random-2-test
  (testing "x^4 + 2x^3 - 15x^2
            results = -5, 0, 3"
    (let [f {4 1.0,
             3 2.0,
             2 -15.0}
          result #{-5.0 0.0 3.0}]
      (is (= result (newton-roots f))))))

(deftest lesson-1-test
  (testing "x^3 - 1
            results = 1"
    (let [f {3 1.0,
             0 -1.0}
          result #{1.0}]
      (is (= result (newton-roots f))))))

(deftest lesson-2-test
  (testing "x^3 + 2
            results = -1.26"
    (let [f {3 1.0,
             0 2.0}
          result #{-1.2599}]
      (is (= result (newton-roots f))))))

(deftest lesson-3-test
  (testing "x^3 - 2.3x^2 + 1.32x
            results = 1.1, 1.2, 0"
    (let [f {3 1.0,
             2 -2.3,
             1 1.32}
          result #{1.1 1.2 0.0}]
      (is (= result (newton-roots f))))))

(deftest lesson-4-test
  (testing "x^3 + 520x^2 - 0.25x - 130
            results = -520, -0.5, 0.5"
    (let [f {3 1.0,
             2 520.0,
             1 -0.25,
             0 -130.0}
          result #{-520.0 -0.5 0.5}]
      (is (= result (newton-roots f))))))

(deftest lesson-5-test
  (testing "25x^3 - 170x^2 - 736x - 640
            results = -1.6, 10"
    (let [f {3 25.0,
             2 -170.0,
             1 -736.0,
             0 -640.0}
          result #{-1.6 10.0}]
      (is (= result (newton-roots f))))))

(deftest lesson-7-test
  (testing "x^4 - 9.2x^3 + 30.2x^2 - 41.2x + 19.2
            results = 1, 2, 3, 3.2"
    (let [f {4 1.0,
             3 -9.2,
             2 30.2,
             1 -41.2,
             0 19.2}
          result #{1.0 2.0 3.0 3.2}]
      (is (= result (newton-roots f))))))

(deftest lesson-8-test
  (testing "x^4 - 6x^3 - 11x^2 + 60x + 100
            results = -2, 5"
    (let [f {4 1.0,
             3 -6.0,
             2 -11.0,
             1 60.0,
             0 100.0}
          result #{-2.0 5.0}]
      (is (= result (newton-roots f))))))
