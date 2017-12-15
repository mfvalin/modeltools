VERSION = 1.0

EXEC = run_work_stream_$(VERSION) SubmitJobContainers_$(VERSION) pick_and_run_work_$(VERSION) r.remove_barrier r.barrier pbs-ompi-sjob-cfg.sh

DOTTED = r.makebarrier.dot

DEMO = multi_container.cfg who-where.sh job_modela  job_modelb  job_post who-where.sh the_jobs_multi

PACKAGE = pbs-ompi-sjob_$(VERSION)_all

TARGETS = $(EXEC) $(DOTTED)

ssm: $(TARGETS) $(DEMO) $(PACKAGE)
	cp $(TARGETS) $(PACKAGE)/bin
	cd $(PACKAGE)/bin ; chmod 755 $(EXEC) ; chmod 644 $(DOTTED)
	mkdir -p $(PACKAGE)/share/example
	cp $(DEMO) $(PACKAGE)/share/example
	tar zcf $(PACKAGE).ssm $(PACKAGE)

dist:	ssm
	cp $(PACKAGE).ssm ~/repository/.

$(PACKAGE):
	s.ssm-prep -D $(PACKAGE) 2>/dev/null

clean:
	rm -rf $(PACKAGE)

distclean: clean
	rm -f $(PACKAGE).ssm
