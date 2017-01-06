(ns ^{:author "Kevin Stadler"}
  usm.measures
  (:require
    [incanter.stats :refer [mean]]))

; population measures

(defn defaultmeasures [x]
  [(mean x) (apply min x) (apply max x)])
(def defaultmeasurenames ["" "min" "max"])

; default bindings
(def ^:dynamic measurefn defaultmeasures)
(def ^:dynamic measurenames defaultmeasurenames)
(def ^:dynamic observedvars ["x" "m" "int"])

; new
(def ^:dynamic xmeasures defaultmeasures)
(def ^:dynamic mmeasures defaultmeasures)
(def ^:dynamic intmeasures defaultmeasures)

(defmacro with-discrete-momentum-measures [& body]
  `(let [orig# measurefn]
     (with-redefs [measurefn #(conj (orig# %) (mean (map (fn [x#] (Math/signum (double x#))) %)))
                   measurenames (conj measurenames "sign")]
       ~@body)))


; simulation run summary measures

(defn has-transition [data]
  (and
    (some #(> 0.5 %) data)
    (some #(< 0.5 %) data)))

(defn indexof
  "Returns the index (starting from 0) of the first item in the collection for
   which (pred item) returns true"
  [pred coll]
  (first (keep-indexed #(when (pred %2) %1) coll)))
;(indexof (partial <= 3) (range 10))
