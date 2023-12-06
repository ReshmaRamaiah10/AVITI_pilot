HTSLIB = htslib-1.17
HTSHEADERS = $(HTSLIB)/htslib/bgzf.h $(HTSLIB)/htslib/hts.h $(HTSLIB)/htslib/sam.h $(HTSLIB)/htslib/thread_pool.h
LIBHTS = $(HTSLIB)/libhts.a

# Compilation variables
CXX = g++
CXXFLAGS = -std=gnu++11 -Wall -I. -I$(HTSLIB) -O3

OBJS1 = RefSeq.o Refs.o demo.o
PROGS1 = demo

PROGRAMS = $(PROGS1) 

.PHONY : all  clean

all : $(PROGRAMS)


$(LIBHTS) :
	cd $(HTSLIB) && $(MAKE)

# Compile objects
$(OBJS1) :
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Generate executables
$(PROGS1) :
	$(CXX) -o $@ $^ -lz -lm -lbz2 -llzma -lcurl -lpthread



# Clean
clean :
	rm -f *.o *~ $(PROGRAMS)
	cd $(HTSLIB) && $(MAKE) clean

# Dependencies for executables
demo : demo.o RefSeq.o Refs.o $(LIBHTS) SEQstring.o SamHeaderText.o SamParser.o BamWriter.o BamAlignment.o

# Dependencies for objects
RefSeq.o : RefSeq.cpp utils.h my_assert.h RefSeq.hpp
Refs.o : Refs.cpp utils.h my_assert.h RefSeq.hpp Refs.hpp
SEQstring.o : SEQstring.cpp $(HTSHEADERS) SEQstring.hpp
SamHeaderText.o : SamHeaderText.cpp $(HTSHEADERS) my_assert.h SamHeaderText.hpp
SamParser.o : SamParser.cpp $(HTSHEADERS) my_assert.h SamHeaderText.hpp SamParser.hpp
BamWriter.o : BamWriter.cpp $(HTSHEADERS) my_assert.h SamHeaderText.hpp BamWriter.hpp
BamAlignment.o : BamAlignment.cpp $(HTSHEADERS) utils.h my_assert.h CIGARstring.hpp SEQstring.hpp QUALstring.hpp MDstring.hpp SamHeaderText.hpp SamParser.hpp BamWriter.hpp BamAlignment.hpp
demo.o : demo.cpp utils.h my_assert.h RefSeq.hpp Refs.hpp $(HTSHEADERS) utils.h my_assert.h CIGARstring.hpp SEQstring.hpp QUALstring.hpp MDstring.hpp SamHeaderText.hpp SamParser.hpp BamWriter.hpp BamAlignment.hpp





# # PROGS2 = rsem-simulate-reads rsem-parse-alignments rsem-run-em rsem-tbam2gbam rsem-bam2wig rsem-bam2readdepth rsem-get-unique rsem-sam-validator rsem-scan-for-paired-end-reads
# # PROGS3 = rsem-run-gibbs rsem-calculate-credibility-intervals










# # $(OBJS2) :
# # 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -O2 -c -o $@ $<

# $(OBJS2) :
# 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -O3 -ffast-math -c -o $@ $<



# $(PROGS2) :
# 	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) -lbz2 -llzma -lcurl -lcrypto -lm -lz -lpthread

# $(PROGS3) :
# 	$(CXX) $(LDFLAGS) -pthread -o $@ $^ $(LDLIBS)


# rsem-gbam2tbam : gbam2tbam.o $(SAMLIBS) Transcript.o Transcripts.o GenomeMap.o SEQstring.o SamHeaderText.o SamParser.o BamWriter.o BamAlignment.o
# rsem-tbam2gbam : tbam2gbam.o $(SAMLIBS) Transcript.o Transcripts.o SEQstring.o SamHeaderText.o SamParser.o BamWriter.o BamAlignment.o



# gbam2tbam.o : gbam2tbam.cpp $(SAMHEADERS) utils.h my_assert.h Transcript.hpp Transcripts.hpp GenomeMap.hpp CIGARstring.hpp SEQstring.hpp QUALstring.hpp MDstring.hpp SamHeaderText.hpp SamParser.hpp BamWriter.hpp BamAlignment.hpp AlignmentGroup.hpp ConversionGroup.hpp
# tbam2gbam.o : tbam2gbam.cpp $(SAMHEADERS) utils.h my_assert.h Transcript.hpp Transcripts.hpp CIGARstring.hpp SEQstring.hpp QUALstring.hpp MDstring.hpp SamHeaderText.hpp SamParser.hpp BamWriter.hpp BamAlignment.hpp AlignmentGroup.hpp ConversionGroup.hpp
# GroupInfo.o : GroupInfo.cpp my_assert.h GroupInfo.hpp
# Orientation.o : Orientation.cpp utils.h sampling.hpp Orientation.hpp
# NoiseProfile.o : NoiseProfile.cpp utils.h $(BOOST)/random.hpp sampling.hpp SEQstring.hpp NoiseProfile.hpp
# NoiseQProfile.o : NoiseQProfile.cpp utils.h $(BOOST)/random.hpp sampling.hpp SEQstring.hpp QUALstring.hpp NoiseProfile.hpp
# Profile.o : Profile.cpp utils.h $(BOOST)/random.hpp sampling.hpp Profile.hpp
# QProfile.o : QProfile.cpp utils.h $(BOOST)/random.hpp sampling.hpp QProfile.hpp
# QualDist.o : QualDist.cpp utils.h $(BOOST)/random.hpp sampling.hpp QUALstring.hpp QualDist.hpp
# Markov.o : Markov.cpp utils.h $(BOOST)/random.hpp sampling.hpp Markov.hpp
# MateLenDist.o : MateLenDist.cpp utils.h $(BOOST)/random.hpp sampling.hpp MateLenDist.hpp
# IlluminaSequenceModel.o: IlluminaSequenceModel.cpp utils.h my_assert.h $(BOOST)/random.hpp sampling.hpp CIGARstring.hpp MDstring.hpp SEQstring.hpp QUALstring.hpp MateLenDist.hpp Markov.hpp Profile.hpp QProfile.hpp QualDist.hpp NoiseProfile.hpp NoiseQProfile.hpp IlluminaSequenceModel.hpp
# FragLenDist.o :

# # wiggle.o: wiggle.cpp $(SAMHEADERS) sam_utils.h utils.h my_assert.h wiggle.h
# # bam2wig.o : bam2wig.cpp utils.h my_assert.h wiggle.h
# # bam2readdepth.o : bam2readdepth.cpp utils.h my_assert.h wiggle.h
# # getUnique.o : getUnique.cpp $(SAMHEADERS) sam_utils.h utils.h 
# # samValidator.o : samValidator.cpp $(SAMHEADERS) sam_utils.h utils.h my_assert.h
# # scanForPairedEndReads.o : scanForPairedEndReads.cpp $(SAMHEADERS) sam_utils.h utils.h my_assert.h 

# # EM.o : EM.cpp $(SAMHEADERS) utils.h my_assert.h Read.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h SingleHit.h PairedEndHit.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h Refs.h GroupInfo.h HitContainer.h ReadIndex.h ReadReader.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h ModelParams.h RefSeq.h RefSeqPolicy.h PolyARules.h Profile.h NoiseProfile.h Transcript.h Transcripts.h HitWrapper.h BamWriter.h simul.h sam_utils.h SamHeader.hpp sampling.h $(BOOST)/boost/random.hpp WriteResults.h
# # Gibbs.o : Gibbs.cpp utils.h my_assert.h $(BOOST)/boost/random.hpp sampling.h simul.h Read.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h SingleHit.h PairedEndHit.h ReadIndex.h ReadReader.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h Profile.h NoiseProfile.h ModelParams.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h RefSeq.h RefSeqPolicy.h PolyARules.h Refs.h GroupInfo.h WriteResults.h 
# # calcCI.o : calcCI.cpp utils.h my_assert.h $(BOOST)/boost/random.hpp sampling.h simul.h Read.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h SingleHit.h PairedEndHit.h ReadIndex.h ReadReader.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h Profile.h NoiseProfile.h ModelParams.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h RefSeq.h RefSeqPolicy.h PolyARules.h Refs.h GroupInfo.h WriteResults.h Buffer.h 
# # simulation.o : simulation.cpp utils.h Read.h SingleRead.h SingleReadQ.h PairedEndRead.h PairedEndReadQ.h Model.h SingleModel.h SingleQModel.h PairedEndModel.h PairedEndQModel.h Refs.h RefSeq.h GroupInfo.h Transcript.h Transcripts.h Orientation.h LenDist.h RSPD.h QualDist.h QProfile.h NoiseQProfile.h Profile.h NoiseProfile.h simul.h $(BOOST)/boost/random.hpp WriteResults.h

# # Dependencies for header files, only include direct headers
# Transcripts.hpp : Transcript.hpp
# RefSeq.hpp : utils.h my_assert.h CIGARstring.hpp MDstring.hpp SEQstring.hpp
# Refs.hpp : RefSeq.hpp
# SamHeaderText.hpp: $(SAMHEADERS) my_assert.h
# SamParser.hpp : $(SAMHEADERS) SamHeaderText.hpp
# BamWriter.hpp : $(SAMHEADERS) my_assert.h SamHeaderText.hpp
# BamAlignment.hpp : $(SAMHEADERS) my_assert.h CIGARstring.hpp SEQstring.hpp QUALstring.hpp MDstring.hpp SamParser.hpp BamWriter.hpp
# AlignmentGroup.hpp : $(SAMHEADERS) SEQstring.hpp QUALstring.hpp SamParser.hpp BamWriter.hpp BamAlignment.hpp
# ConversionGroup.hpp : BamAlignment.hpp AlignmentGroup.hpp
# sampling.hpp : $(BOOST)/random.hpp
# Orientation.hpp : utils.h sampling.hpp
# NoiseProfile.hpp: utils.h sampling.hpp SEQstring.hpp
# NoiseQProfile.hpp: utils.h sampling.hpp SEQstring.hpp QUALstring.hpp
# Profile.hpp: utils.h sampling.hpp
# QProfile.hpp: utils.h sampling.hpp
# QualDist.hpp: utils.h sampling.hpp QUALstring.hpp
# Markov.hpp: utils.h sampling.hpp
# MateLenDist.hpp : utils.h sampling.hpp
# IlluminaSequenceModel.hpp : utils.h sampling.hpp RefSeq.hpp CIGARstring.hpp QUALstring.hpp SEQstring.hpp MateLenDist.hpp Markov.hpp Profile.hpp QProfile.hpp QualDist.hpp NoiseProfile.hpp NoiseQProfile.hpp
# FragLenDist.hpp : 

