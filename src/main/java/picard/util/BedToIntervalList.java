package picard.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.IntervalsManipulationProgramGroup;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * @author nhomer
 */
@CommandLineProgramProperties(
        summary = BedToIntervalList.USAGE_SUMMARY + BedToIntervalList.USAGE_DETAILS,
        oneLineSummary = BedToIntervalList.USAGE_SUMMARY,
        programGroup = IntervalsManipulationProgramGroup.class
)
@DocumentedFeature
public class BedToIntervalList extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Converts a BED file to a Picard Interval List.  " ;
    static final String USAGE_DETAILS = "This tool provides easy conversion from BED to the Picard interval_list format which is " +
            "required by many Picard processing tools. Note that the coordinate system of BED files is such that the first base or " +
            "position in a sequence is numbered \"0\", while in interval_list files it is numbered \"1\"." +
            "<br /><br />" +
            "BED files contain sequence data displayed in a flexible format that includes nine optional fields, " +
            "in addition to three required fields within the annotation tracks. The required fields of a BED file include:" +
            "<pre>" +
            "     chrom - The name of the chromosome (e.g. chr20) or scaffold (e.g. scaffold10671) <br />"   +
            "     chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered \"0\" <br />"   +
            "     chromEnd - The ending position of the feature in the chromosome or scaffold.  The chromEnd base is not" +
            " included in the display of the feature. For example, the first 100 bases of a " +
            "chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99." +
            "</pre>" +
            "In each annotation track, the number of fields per line must be consistent throughout a data set. " +
            "For additional information regarding BED files and the annotation field options, please see:" +
            " http://genome.ucsc.edu/FAQ/FAQformat.html#format1." +
            "<br /> <br /> " +
            "Interval_list files contain sequence data distributed into intervals. The interval_list file format is relatively simple " +
            "and reflects the SAM alignment format to a degree.  A SAM style header must be present in the file that lists the sequence " +
            "records against which the intervals are described.  After the header, the file then contains records, one per line in plain " +
            "text format with the following values tab-separated::" +
            "<pre> " +
            "     -Sequence name (SN) - The name of the sequence in the file for identification purposes, can be chromosome number e.g. chr20 <br /> " +
            "     -Start position - Interval start position (starts at +1) <br /> " +
            "     -End position - Interval end position (1-based, end inclusive) <br /> " +
            "     -Strand - Indicates +/- strand for the interval (either + or -) <br /> " +
            "     -Interval name - (Each interval should have a unique name) " +
            "</pre>" +
            "<br/>" +
            "This tool requires a sequence dictionary, provided with the SEQUENCE_DICTIONARY or SD argument. " +
            "The value given to this argument can be any of the following:" +
            "<pre>" +
            "    - A file with .dict extension generated using Picard's CreateSequenceDictionaryTool</br>" +
            "    - A reference.fa or reference.fasta file with a reference.dict in the same directory</br>" +
            "    - Another IntervalList with @SQ lines in the header from which to generate a dictionary</br>" +
            "    - A VCF that contains #contig lines from which to generate a sequence dictionary</br>" +
            "    - A SAM or BAM file with @SQ lines in the header from which to generate a dictionary</br>" +
            "</pre>" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar BedToIntervalList \\<br />" +
            "      I=input.bed \\<br />" +
            "      O=list.interval_list \\<br />" +
            "      SD=reference_sequence.dict" +
            "</pre>" +
            "<br /> <br /> "+
            "<hr />"
            ;
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input BED file")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME,
            doc = "The sequence dictionary, or BAM/VCF/IntervalList from which a dictionary can be extracted.")
    public File SEQUENCE_DICTIONARY;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output Picard Interval List")
    public File OUTPUT;

    @Argument(doc="If true, sort the output interval list before writing it.")
    public boolean SORT = true;

    @Argument(doc="If true, unique the output interval list by merging overlapping regions, before writing it (implies sort=true).")
    public boolean UNIQUE = false;

    @Hidden
    @Argument(doc = "If true, entries that are on contig-names that are missing from the provided dictionary will be dropped.")
    public boolean DROP_MISSING_CONTIGS = false;

    @Argument(doc = "If true, write length zero intervals in input bed file to resulting interval list file.")
    public boolean KEEP_LENGTH_ZERO_INTERVALS = false;

    private final Log LOG = Log.getInstance(getClass());

    @Override
    protected int doWork() {
        // Only assert readability for regular files; FIFOs, named pipes, /dev/stdin,
        // and process-substitution paths (/dev/fd/N) all return false from isFile().
        if (INPUT.isFile()) {
            IOUtil.assertFileIsReadable(INPUT);
        }
        IOUtil.assertFileIsReadable(SEQUENCE_DICTIONARY);
        IOUtil.assertFileIsWritable(OUTPUT);

        try {
            final SAMFileHeader header = new SAMFileHeader();
            final SAMSequenceDictionary samSequenceDictionary = SAMSequenceDictionaryExtractor.extractDictionary(SEQUENCE_DICTIONARY.toPath());
            header.setSequenceDictionary(samSequenceDictionary);
            header.setSortOrder(SAMFileHeader.SortOrder.coordinate);

            // For /dev/stdin specifically, wrap System.in so that tests can inject data
            // via System.setIn().  All other non-regular paths (named pipes, /dev/fd/N,
            // process substitutions) are opened by path through FileReader as normal.
            try (final BufferedReader reader = INPUT.getPath().equals("/dev/stdin")
                    ? new BufferedReader(new InputStreamReader(System.in))
                    : new BufferedReader(new FileReader(INPUT))) {
                // Sniff the format before parsing; reject anything that isn't BED.
                reader.mark(8 * 1024);
                final FormatDetectionResult detected = detectIntervalFormat(reader);
                if (detected.format() != IntervalFileFormat.BED) {
                    final String hint = detected.format() == IntervalFileFormat.INTERVAL_LIST
                            ? " Input appears to be an interval_list file; supply a BED file instead."
                            : (detected.firstLine() != null ? " First data line: " + detected.firstLine() : " File appears to be empty or contain only headers.");
                    throw new PicardException("BedToIntervalList requires BED format input." + hint);
                }
                reader.reset();
                IntervalList out = fromBed(reader, header, DROP_MISSING_CONTIGS, KEEP_LENGTH_ZERO_INTERVALS, LOG);
                if (SORT) out = out.sorted();
                if (UNIQUE) out = out.uniqued();
                out.write(OUTPUT);
                LOG.info(String.format("Wrote %d intervals spanning a total of %d bases",
                        out.getIntervals().size(), out.getBaseCount()));
            }
        } catch (final IOException e) {
            throw new RuntimeException(e);
        }

        return 0;
    }

    // -------------------------------------------------------------------------
    // Format-sniffing / generic interval loading
    // -------------------------------------------------------------------------

    private enum IntervalFileFormat { INTERVAL_LIST, BED, UNKNOWN }

    private record FormatDetectionResult(IntervalFileFormat format, String firstLine) {}

    /**
     * Detects whether a reader contains interval_list or BED content by inspecting the first
     * significant (non-empty, non-comment) line.  The reader's position is NOT reset after this
     * call — callers must {@link BufferedReader#mark} and {@link BufferedReader#reset} as needed.
     *
     * <ul>
     *   <li>Lines starting with {@code @} indicate interval_list (SAM-style header).</li>
     *   <li>Lines with ≥3 tab-separated fields indicate BED.</li>
     *   <li>Comment lines starting with {@code #} are skipped.</li>
     * </ul>
     */
    private static FormatDetectionResult detectIntervalFormat(final BufferedReader reader) throws IOException {
        String line;
        while ((line = reader.readLine()) != null) {
            final String trimmed = line.trim();
            if (trimmed.isEmpty() || trimmed.startsWith("#")) {
                continue;
            }
            if (trimmed.startsWith("@")) {
                return new FormatDetectionResult(IntervalFileFormat.INTERVAL_LIST, null);
            } else if (trimmed.split("\t").length >= 3) {
                return new FormatDetectionResult(IntervalFileFormat.BED, null);
            } else {
                return new FormatDetectionResult(IntervalFileFormat.UNKNOWN, trimmed);
            }
        }
        return new FormatDetectionResult(IntervalFileFormat.UNKNOWN, null);
    }

    /**
     * Loads an {@link IntervalList} from {@code reader}, auto-detecting whether the content is
     * BED or interval_list format.  The reader must support {@link BufferedReader#mark} /
     * {@link BufferedReader#reset} (a plain {@link BufferedReader} always does — this works for
     * regular files, pipes, FIFOs, and process-substitution streams alike).
     *
     * <p>Returns a sorted, uniqued {@link IntervalList}.
     *
     * @param reader     source; must support mark/reset
     * @param dictionary sequence dictionary used when parsing BED-format content
     * @throws IOException  on read error
     * @throws PicardException on unrecognized format
     */
    public static IntervalList loadIntervals(final BufferedReader reader,
                                             final SAMSequenceDictionary dictionary) throws IOException {
        // 8 KB is enough to cover any realistic BED/interval_list preamble.
        // BufferedReader.mark() is backed by an in-memory char array, so this also works
        // correctly for non-seekable sources such as pipes and FIFOs.
        reader.mark(8 * 1024);
        final FormatDetectionResult detected = detectIntervalFormat(reader);
        reader.reset();
        return switch (detected.format()) {
            case INTERVAL_LIST -> IntervalList.fromReader(reader).uniqued();
            case BED -> {
                final SAMFileHeader header = new SAMFileHeader();
                header.setSequenceDictionary(dictionary);
                yield fromBed(reader, header, false, false, null).uniqued();
            }
            case UNKNOWN -> throw new PicardException(
                    "Unrecognized interval file format. Expected interval_list (lines starting with @) " +
                    "or BED (≥3 tab-separated fields). First data line: " + detected.firstLine());
        };
    }

    /**
     * Convenience overload of {@link #loadIntervals(BufferedReader, SAMSequenceDictionary)} that
     * opens {@code file} itself.  {@link IOException} is wrapped in a {@link PicardException}.
     *
     * @param file       BED or interval_list file to read
     * @param dictionary sequence dictionary used when parsing BED-format content
     */
    public static IntervalList loadIntervals(final File file, final SAMSequenceDictionary dictionary) {
        try (final BufferedReader reader = new BufferedReader(new FileReader(file))) {
            return loadIntervals(reader, dictionary);
        } catch (final IOException e) {
            throw new PicardException("Error reading intervals from " + file, e);
        }
    }

    // -------------------------------------------------------------------------
    // BED-only parsing
    // -------------------------------------------------------------------------

    /**
     * Parses BED records from a {@link BufferedReader} into an {@link IntervalList}.
     * Comment lines starting with {@code #} are skipped.
     * Coordinates are converted from 0-based half-open BED to 1-based closed intervals.
     * <p>
     * This method works with any reader, including streams backed by pipes or FIFOs.
     *
     * @param reader                  source of BED records
     * @param header                  SAMFileHeader whose sequence dictionary is used for validation
     * @param dropMissingContigs      if true, records on contigs absent from the dictionary are silently
     *                                dropped; if false, such records throw a {@link PicardException}
     * @param keepLengthZeroIntervals if true, length-zero BED intervals (chromStart == chromEnd) are
     *                                included; if false, they are silently dropped
     * @param log                     used for informational and warning messages; may be {@code null}
     *                                to suppress all logging
     * @return an unsorted, non-uniqued {@link IntervalList}
     * @throws IOException on read error
     */
    public static IntervalList fromBed(final BufferedReader reader,
                                       final SAMFileHeader header,
                                       final boolean dropMissingContigs,
                                       final boolean keepLengthZeroIntervals,
                                       final Log log) throws IOException {
        final IntervalList intervalList = new IntervalList(header);
        int missingIntervals = 0;
        int missingRegion = 0;
        int lengthZeroIntervals = 0;

        String line;
        while ((line = reader.readLine()) != null) {
            final String trimmed = line.trim();
            if (trimmed.isEmpty() || trimmed.startsWith("#")) {
                continue;
            }

            final String[] fields = trimmed.split("\t");
            if (fields.length < 3) {
                throw new PicardException("Invalid BED line (fewer than 3 tab-separated fields): " + line);
            }

            final String sequenceName = fields[0];
            final int start = Integer.parseInt(fields[1]) + 1; // BED start is 0-based; convert to 1-based
            final int end   = Integer.parseInt(fields[2]);     // BED end is 0-based exclusive == 1-based inclusive
            // NB: do not use an empty name within an interval
            final String name = (fields.length > 3 && !fields[3].isEmpty()) ? fields[3] : null;
            final boolean isNegativeStrand = fields.length > 5 && "-".equals(fields[5]);

            final SAMSequenceRecord sequenceRecord = header.getSequenceDictionary().getSequence(sequenceName);

            if (null == sequenceRecord) {
                if (dropMissingContigs) {
                    if (log != null) log.info(String.format("Dropping interval with missing contig: %s:%d-%d", sequenceName, start, end));
                    missingIntervals++;
                    missingRegion += end - start + 1;
                    continue;
                }
                throw new PicardException(String.format("Sequence '%s' was not found in the sequence dictionary", sequenceName));
            } else if (start < 1) {
                throw new PicardException(String.format("Start on sequence '%s' was less than one: %d", sequenceName, start));
            } else if (sequenceRecord.getSequenceLength() < start) {
                throw new PicardException(String.format("Start on sequence '%s' was past the end: %d < %d", sequenceName, sequenceRecord.getSequenceLength(), start));
            } else if ((end == 0 && start != 1) || end < 0) {
                throw new PicardException(String.format("End on sequence '%s' was less than one: %d", sequenceName, end));
            } else if (sequenceRecord.getSequenceLength() < end) {
                throw new PicardException(String.format("End on sequence '%s' was past the end: %d < %d", sequenceName, sequenceRecord.getSequenceLength(), end));
            } else if (end < start - 1) {
                throw new PicardException(String.format("On sequence '%s', end < start-1: %d <= %d", sequenceName, end, start));
            }

            if ((start == end + 1) && !keepLengthZeroIntervals) {
                if (log != null) log.info(String.format("Skipping writing length zero interval at %s:%d-%d.", sequenceName, start, end));
                lengthZeroIntervals++;
                continue;
            }
            if (start == end + 1) {
                lengthZeroIntervals++;
            }

            intervalList.add(new Interval(sequenceName, start, end, isNegativeStrand, name));
        }

        if (log != null) {
            if (dropMissingContigs) {
                if (missingRegion == 0) {
                    log.info("There were no missing regions.");
                } else {
                    log.warn(String.format("There were %d missing regions with a total of %d bases", missingIntervals, missingRegion));
                }
            }
            if (!keepLengthZeroIntervals) {
                if (lengthZeroIntervals == 0) {
                    log.info("No input regions had length zero, so none were skipped.");
                } else {
                    log.info(String.format("Skipped writing a total of %d entries with length zero in the input file.", lengthZeroIntervals));
                }
            } else {
                if (lengthZeroIntervals > 0) {
                    log.warn(String.format("Input file had %d entries with length zero. Run with the KEEP_LENGTH_ZERO_INTERVALS flag set to false to remove these.", lengthZeroIntervals));
                }
            }
        }

        return intervalList;
    }
}
