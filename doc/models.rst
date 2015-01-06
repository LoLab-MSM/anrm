Overview of the models
======================

The Apoptosis Necroptosis Reaction Network (ANRM) contains roughly four
major mechanistic groups. Each mechanistic group employs modules that list the
reaction steps  (or rules) contained group.

- TNF ligation and Complex I assembly
- Complex II reactions
- Hypothesized Bid regulation of RIP1 activity
- Mitochondrial Apoptotic Signaling and PARP Cleavage.

.. image:: images/ANRM_pathway1.png
   :height: 400
   :align: center
   :alt: EARM 2.0 pathway/architecture diagram

The apoptosis-necroptosis decision-making capability of the ANRM models is captured 
primarily in the assembly and activities of Complex II. There exists distinct RIP1 
and TRADD dependent complexes IIa and IIb respectively. Apoptosis is potentiated by 
RIP1 dependent Complex IIa, which catalyzes the activation of pro-apoptotic caspase 
8. When caspase 8 activity is suppressed (e.g. by cFlip-S) phosphorylated RIP1 
(RIP1-po4) molecules can accumulate trigger activation of pro-necroptotic MLKL. 
TRADD dependent complex IIb mediates the RIP1 dependent FADD independent 
necroptosis observed by [Cite2014]_. Complex IIb can, in our model, however, 
recruit FADD and concomitantly signal apoptosis.

The initial and primary purpose of the ANRM is to explore hypothesized roles of Bid
in regulating and RIP1 mediated necroptosis. Zinkel et al. proposed that 
phosphorylated Bid can recruit and inactivate (via truncation) RIP1.  This model is 
described in the **"irvin_anrm_alt_hypothesis"** model. We refined this model to
include a complex that contains phosphorylated Bid and the Complex II dependent 
procaspase 8/cFlip-L heterodimer. This complex  catalyzed RIP1 truncation. The 
refined model is model described in **"irvin_anrm_model"**. A model lacking
phosphorylated Bid is described in **"irvin_anrm_wo_bidpo4"**.

Additional hypotheses were proposed in order to address a knowledge gap in the 
transition from Complex I on the plasma membrane to the cytoplasmic Complex II. 
Current descriptions state that complex I recruits, poly-ubiquitinates, de-
ubiquitinates, and releases RIP1 into the cytoplasm where it will trigger 
programmed cell death. These descriptions describe a net-zero change in the state 
of RIP1. We addressed the problem by hypothesizing various modifications that might 
distinguish endogenous RIP1 from that released from Complex I. Three models 
accurately reproduce the qualitative behavior of apoptotic/necroptotic cells:

- FADD transiently localizes to TNFR1 to retrieve RIP1
- TRADD:RIP1 is released from Complex I and FADD replaces TRADD in TRADD:RIP1.
- TRADD:RIP1 is released from Complex I and FADD binds TRADD in TRADD:RIP1.

These three models can be superimposed for a total of 7 possible mechanisms. The 
ANRM models use a superpostion of all three. Combining the Bid hypotheses and the 
models addressing the transition from Complex I to Complex II results in 21 total 
models.

Mitochondrial Apoptotic Signaling and PARP Cleavage pathway components and reaction 
topologies  are re-used from the previously published EARM 2.0 [Albeck2010]_.

The models in EARM
------------------

Below is a list of the 15 alternative Bcl-2 reaction topologies incorporated
into EARM. More detailed descriptions of each model, along with the source
code, are found in :doc:`modules/index`.

- M1a/b: EARM 2.0, Embedded [Lopez2013]_ 
- M2a/b: EARM 2.0, Direct [Lopez2013]_
- M3a/b: EARM 2.0, Indirect [Lopez2013]_
- M4a/b: "Minimal Model" (Figure 11b) from Albeck et al. (2008) [Albeck2008]_
- M5a/b: "Model B + Bax multimerization" (Figure 11c) from Albeck et al. (2008)
  [Albeck2008]_
- M6a/b: "Model C + mitochondrial transport" (Figure 11d) from Albeck et al.
  (2008) [Albeck2008]_
- M7a/b: "Current model" (Figure 11e) from Albeck et al. (2008) [Albeck2008]_
- M8a/b: "Current model + cooperativity" (Figure 11f) from Albeck et al. (2008)
  [Albeck2008]_
- M9a/b: Deterministic model from Chen et al. (2007), Biophysical Journal
  [Chen2007biophysj]_
- M10a/b: Indirect model from Chen et al. (2007), FEBS Letters [Chen2007febs]_
- M11a/b: Direct model from Chen et al. (2007), FEBS Letters [Chen2007febs]_
- M12a/b: Direct model from Cui et al. (2008) [Cui2008]_
- M13a/b: Direct model 1 from Cui et al. (2008) [Cui2008]_
- M14a/b: Direct model 2 from Cui et al. (2008) [Cui2008]_
- M15a/b: Model incorporating Bad phosphorylation from Howells et al. (2011)
  [Howells2011]_

How the model code is organized
-------------------------------

For each of the 30 models, there is a corresponding ``.py`` file containing
the model definition. This way any model can be imported using the
straightforward syntax (for example, for model M1a)::

    from earm.lopez_embedded import model

However, the Python files for each individual model in general do not contain
much code--they mainly call functions from other modules. For example,
here is the source code for the file ``earm/lopez_embedded.py``, which
implements model M1a:

.. literalinclude:: ../earm/lopez_embedded.py

As this example shows, the model file calls a series of macros that declare the
monomers (:py:func:`earm.albeck_modules.ligand_to_c8_monomers`,
:py:func:`earm.lopez_modules.momp_monomers`, and
:py:func:`earm.albeck_modules.apaf1_to_parp_monomers`), then calls the macros
implementing the upstream and downstream pathway elements
(:py:func:`earm.albeck_modules.rec_to_bid` and
:py:func:`earm.albeck_modules.pore_to_parp`), and finally calls the macro for
the specific Bcl-2 topology involved: :py:func:`earm.lopez_modules.embedded`.
Since the observables for all of the full apoptosis model variants is the same,
these are declared in the final macro that is called,
:py:func:`earm.shared.observables`.

All of the model ``.py`` files follow this pattern, calling a handful of macros
to declare monomers and observables and select implementations of different
pathway modules.

Note that the ``.py`` model files for the full-apoptosis models (M1a - M15a) 
are found in the top-level ``earm`` module, but the files for the MOMP-only
models (M1b - M15b) are found in the submodule ``earm.mito``.

The documentation for all model files (with links to source code) can be found
at the following links:

.. toctree::
   :maxdepth: 2

   modules/model_list

The details of the various macro implementations are found in the following
four files:

.. toctree::
   :maxdepth: 2

   modules/lopez_modules
   modules/albeck_modules
   modules/shen_modules
   modules/shared

MOMP module "boundaries"
------------------------

In the interest of consistency, all of the models have been defined with the
same boundaries in terms of their position in the overall extrinsic apoptosis
pathway: they are all triggered by the addition of an active BH3-only species
(e.g., tBid) as their most "upstream" event, and they all result in the release
in one or more mitochondrial substances (e.g. Cytochrome C and/or Smac) as
their most downstream event. This represents a compromise between the approach
of the MOMP models described in Albeck et al (in which caspase-8, rather than
tBid, served as the input) and the models of the Shen group, in which active
Bax or Bax pores, rather than Cytochrome C or Smac, served as the output.

While these interface boundaries represent the default condition, they can be
modified by passing parameters in to the module macro. For example, by setting
`do_pore_transport=False` in the call to one of the Shen models, the Cytochrome
C and Smac release reactions are not added, and the models can be directly
compared to their originally published versions. Similarly, the upstream
caspase-8/Bid reactions can be added to the Albeck MOMP models to make them
consistent with their published versions.

.. note::MOMP module initial conditions

    The default initial conditions for the MOMP modules is for there to be
    **none of the apoptosis-inducing BH3-only proteins** (i.e., tBid) present
    by default. This means that to reproduce figures from the original
    publications this initial condition will have to be set appropriately.

    If a Bid initial condition is specified, it is for the full-length,
    untruncated form (i.e., Bid(state='U')).

Since our purpose in using these models is primarily to embed them in a common
pathway context, rather than to reproduce previous results for posterity, our
conclusion in working with them was that it is better to have a consistent
interface by default and reproduce published results by modifying the model
rather than implement the model as published by default and then have to
specifically modify each one separately to fit the pathway context
appropriately.

How to use the models
---------------------

To import a model use the syntax::

    from earm.lopez_embedded import model

That's it. You now have a model object that you can query, simulate, perform
parameter estimation on, etc. If you wanted the MOMP-only version, which is 
in the sub-module ``mito``, simply run::

    from earm.mito.lopez_embedded import model

If you want to work with multiple models at the same time (e.g., to compare
them), you can write::

    from earm.chen2007_indirect import model as indirect
    from earm.chen2007_direct import model as direct

For more information on the kinds of analysis you can do using PySB models,
see the `PySB documentation <http://pysb.rtfd.org>`_.

Parameter values
----------------

Parameter values (both rate constants and initial protein concentrations) are
embedded directly in the model code rather than in a separate table or file.
The values in the model definition represent estimates or nominal values and
can be easily overridden using values obtained (for example) by measurement or
parameter estimation algorithms.  We do not maintain a separate list or table of
parameter values, as we have found that the clearest description of the
meaning of a rate parameter is the macro or rule statement in which it is
embedded.

If desired, lists of all model parameters can be obtained via the parameters
instance variable of the model object, i.e.::

    model.parameters

A list of all parameter names can be obtained using the list comprehension::

    [p.name for p in model.parameters]

The code is meant to be read!
-----------------------------

As much as possible, we have attempted to make the code for models themselves
transparent and well-documented. The documentation for each model topology has
been embedded inline in the model source code, and the documentation provided in
the :doc:`modules/index` section of the documentation is drawn directly from
this source.

Moreover, the models have been written using a high-level vocabulary of
frequently re-used macros and motifs, with the aim of revealing broad
similarities and differences between models. The models thus consist of
statements such as::

    translocate_tBid_Bax_BclXL()
    catalyze(Bid(state='T'), Bax(state='M'), Bax(state='A'), klist)

which can be read as saying that "tBid, Bax and BclXL translocate [to the
mitochondrial membrane], and tBid catalyzes Bax from a Mitochondrial (but
inactive) state to an Active state." Understanding the precise mechanisms of
these macros (as expressed in terms of rules and reactions) takes some
familiarity with their implementation, but as there is a fairly limited set
of macros, this should hopefully not present a significant barrier.
