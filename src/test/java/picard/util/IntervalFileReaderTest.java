package picard.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.List;

/**
 * Tests for {@link IntervalFileReader}: format detection, BED parsing, and generic interval loading.
 */
public class IntervalFileReaderTest {

    private static final String TEST_DATA_DIR = "testdata/picard/util/BedToIntervalListTest";

    /** Build a minimal dictionary with chr1..chr8 each of length 1 000 000. */
    private static SAMSequenceDictionary buildDictionary() {
        final SAMSequenceDictionary dict = new SAMSequenceDictionary();
        for (int i = 1; i <= 8; i++) {
            dict.addSequence(new SAMSequenceRecord("chr" + i, 1_000_000));
        }
        return dict;
    }

    private static SAMFileHeader buildHeader() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(buildDictionary());
        return header;
    }

    private static BufferedReader readerOf(final String content) {
        return new BufferedReader(new StringReader(content));
    }

    @Test
    public void testDetectBedFormat() throws IOException {
        try (final BufferedReader reader = readerOf("chr1\t100\t200\n")) {
            final IntervalFileReader.FormatDetectionResult result = IntervalFileReader.detectIntervalFormat(reader);
            Assert.assertEquals(result.format(), IntervalFileReader.IntervalFileFormat.BED);
        }
    }

    @Test
    public void testDetectIntervalListFormat() throws IOException {
        try (final BufferedReader reader = readerOf("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000000\n")) {
            final IntervalFileReader.FormatDetectionResult result = IntervalFileReader.detectIntervalFormat(reader);
            Assert.assertEquals(result.format(), IntervalFileReader.IntervalFileFormat.INTERVAL_LIST);
        }
    }

    @Test
    public void testDetectUnknownFormat() throws IOException {
        try (final BufferedReader reader = readerOf("not_valid_data\n")) {
            final IntervalFileReader.FormatDetectionResult result = IntervalFileReader.detectIntervalFormat(reader);
            Assert.assertEquals(result.format(), IntervalFileReader.IntervalFileFormat.UNKNOWN);
            Assert.assertEquals(result.firstLine(), "not_valid_data");
        }
    }

    @Test
    public void testDetectSkipsComments() throws IOException {
        try (final BufferedReader reader = readerOf("# comment\n\nchr1\t100\t200\n")) {
            final IntervalFileReader.FormatDetectionResult result = IntervalFileReader.detectIntervalFormat(reader);
            Assert.assertEquals(result.format(), IntervalFileReader.IntervalFileFormat.BED);
        }
    }

    @Test
    public void testDetectEmptyFile() throws IOException {
        try (final BufferedReader reader = readerOf("")) {
            final IntervalFileReader.FormatDetectionResult result = IntervalFileReader.detectIntervalFormat(reader);
            Assert.assertEquals(result.format(), IntervalFileReader.IntervalFileFormat.UNKNOWN);
            Assert.assertNull(result.firstLine());
        }
    }

    @Test
    public void testDetectOnlyComments() throws IOException {
        try (final BufferedReader reader = readerOf("# just a comment\n# another\n")) {
            final IntervalFileReader.FormatDetectionResult result = IntervalFileReader.detectIntervalFormat(reader);
            Assert.assertEquals(result.format(), IntervalFileReader.IntervalFileFormat.UNKNOWN);
            Assert.assertNull(result.firstLine());
        }
    }

    @Test
    public void testFromBedSimple() throws IOException {
        final String bed = "chr1\t100\t200\n" +
                           "chr2\t500\t600\n";
        final IntervalList result = IntervalFileReader.fromBed(readerOf(bed), buildHeader(), false, false, null);
        final List<Interval> intervals = result.getIntervals();
        Assert.assertEquals(intervals.size(), 2);
        // BED 0-based half-open -> 1-based closed: (100,200) -> (101,200)
        Assert.assertEquals(intervals.get(0).getContig(), "chr1");
        Assert.assertEquals(intervals.get(0).getStart(), 101);
        Assert.assertEquals(intervals.get(0).getEnd(), 200);
        // (500,600) -> (501,600)
        Assert.assertEquals(intervals.get(1).getContig(), "chr2");
        Assert.assertEquals(intervals.get(1).getStart(), 501);
        Assert.assertEquals(intervals.get(1).getEnd(), 600);
    }

    @Test
    public void testFromBedExtendedFields() throws IOException {
        // 6-field BED with name and strand
        final String bed = "chr1\t100\t2000\tchr1_100_2000+\t11\t+\n" +
                           "chr1\t3000\t4000\tchr1_3000_4000-\t12\t-\n";
        final IntervalList result = IntervalFileReader.fromBed(readerOf(bed), buildHeader(), false, false, null);
        final List<Interval> intervals = result.getIntervals();
        Assert.assertEquals(intervals.size(), 2);
        Assert.assertEquals(intervals.get(0).getName(), "chr1_100_2000+");
        Assert.assertFalse(intervals.get(0).isNegativeStrand());
        Assert.assertEquals(intervals.get(1).getName(), "chr1_3000_4000-");
        Assert.assertTrue(intervals.get(1).isNegativeStrand());
    }

    @Test
    public void testFromBedSkipsCommentsAndBlankLines() throws IOException {
        final String bed = "# comment line\n" +
                           "\n" +
                           "chr1\t100\t200\n" +
                           "# another comment\n" +
                           "chr1\t300\t400\n";
        final IntervalList result = IntervalFileReader.fromBed(readerOf(bed), buildHeader(), false, false, null);
        Assert.assertEquals(result.getIntervals().size(), 2);
    }

    @Test(expectedExceptions = PicardException.class)
    public void testFromBedMissingContigThrows() throws IOException {
        final String bed = "chrX\t100\t200\n";
        IntervalFileReader.fromBed(readerOf(bed), buildHeader(), false, false, null);
    }

    @Test
    public void testFromBedDropMissingContigs() throws IOException {
        final String bed = "chr1\t100\t200\n" +
                           "chrX\t100\t200\n";
        final IntervalList result = IntervalFileReader.fromBed(readerOf(bed), buildHeader(), true, false, null);
        Assert.assertEquals(result.getIntervals().size(), 1);
        Assert.assertEquals(result.getIntervals().get(0).getContig(), "chr1");
    }

    @Test(expectedExceptions = PicardException.class)
    public void testFromBedTooFewFields() throws IOException {
        final String bed = "chr1\t100\n";
        IntervalFileReader.fromBed(readerOf(bed), buildHeader(), false, false, null);
    }

    @Test
    public void testFromBedZeroLengthKept() throws IOException {
        // chromStart == chromEnd means length-zero interval in BED
        final String bed = "chr1\t100\t100\n";
        final IntervalList result = IntervalFileReader.fromBed(readerOf(bed), buildHeader(), false, true, null);
        Assert.assertEquals(result.getIntervals().size(), 1);
        // (100,100) in BED -> start=101, end=100 in 1-based
        Assert.assertEquals(result.getIntervals().get(0).getStart(), 101);
        Assert.assertEquals(result.getIntervals().get(0).getEnd(), 100);
    }

    @Test
    public void testFromBedZeroLengthSkipped() throws IOException {
        final String bed = "chr1\t100\t100\n";
        final IntervalList result = IntervalFileReader.fromBed(readerOf(bed), buildHeader(), false, false, null);
        Assert.assertEquals(result.getIntervals().size(), 0);
    }

    @Test
    public void testLoadIntervalsFromBedReader() throws IOException {
        final String bed = "chr1\t100\t200\n" +
                           "chr1\t300\t400\n";
        final IntervalList result = IntervalFileReader.loadIntervals(readerOf(bed), buildDictionary());
        // loadIntervals returns uniqued result
        Assert.assertEquals(result.getIntervals().size(), 2);
    }

    @Test
    public void testLoadIntervalsFromBedFile() {
        final File bedFile = new File(TEST_DATA_DIR, "simple.bed");
        final IntervalList result = IntervalFileReader.loadIntervals(bedFile, buildDictionary());
        Assert.assertEquals(result.getIntervals().size(), 2);
        // simple.bed: chr1 100 2000  and  chr1 3000 4000
        Assert.assertEquals(result.getIntervals().get(0).getStart(), 101);
        Assert.assertEquals(result.getIntervals().get(0).getEnd(), 2000);
        Assert.assertEquals(result.getIntervals().get(1).getStart(), 3001);
        Assert.assertEquals(result.getIntervals().get(1).getEnd(), 4000);
    }

    @Test
    public void testLoadIntervalsFromIntervalListFile() {
        final File ilFile = new File(TEST_DATA_DIR, "seq_dict_test.dictionary.interval_list");
        // interval_list format doesn't actually need the dictionary arg (it has its own header),
        // but we pass one anyway since loadIntervals requires it for the BED code path.
        final IntervalList result = IntervalFileReader.loadIntervals(ilFile, buildDictionary());
        Assert.assertFalse(result.getIntervals().isEmpty());
    }

    @Test(expectedExceptions = PicardException.class)
    public void testLoadIntervalsUnknownFormatThrows() throws IOException {
        IntervalFileReader.loadIntervals(readerOf("not_valid"), buildDictionary());
    }

    @DataProvider
    public Object[][] outOfBoundsBedData() {
        return new Object[][] {
                // end past contig length
                {"chr1\t0\t1000001\n"},
                // start past contig length
                {"chr1\t1000001\t1000002\n"},
        };
    }

    @Test(dataProvider = "outOfBoundsBedData", expectedExceptions = PicardException.class)
    public void testFromBedOutOfBoundsThrows(final String bed) throws IOException {
        IntervalFileReader.fromBed(readerOf(bed), buildHeader(), false, false, null);
    }
}
